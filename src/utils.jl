# Used to load sequential data where the time index is preceeded by 0's: 0004
function numberDigits( number; maxdigits=4 )
   ndigs = 0
   for d in 0:maxdigits
       ndigs += number ÷ 10^d > 0
   end
   nzeros = ( maxdigits - ndigs );
   return ( nzeros < 0 ) ? string(number) : "0"^(nzeros)*string(number)
end

# computing mean and standard deviation with a single loop 
#
# mean = sum( a )/N
# std = sum( a - mean )^2/N = sum( a*a + mean*mean - 2*a*mean )/N = 
#     	                      sum( a*a )/N + N/N*mean^2 - 2*mean*sum( a )/N =
#                             sum( a*a )/N + mean^2 - 2*mean^2

function meanStd( a::AbstractArray{<:Real,N} ) where {N}

	len   = length(a)
	suma  = 0
	sumaa = 0
	@inbounds @simd for idx in 1:len
		suma  += a[idx]
		sumaa += a[idx]*a[idx]
	end
	mean = suma/len; 
	return mean, sqrt( sumaa/len + mean*mean - 2*mean*mean )
end


function vectorFieldSize( datasize, pp::PIVParameters )

    step = pp.interSize .- pp.overlap
    return length.( StepRange.( pp.interSize, step, datasize ) )
end

# optimized to use SIMD, performing 4 comparisons at the same time
function maxval( a::AbstractArray{<:Real,N} ) where {N}

	len = length(a)

	if len < 4
		return maximum( a ) 
	end

	simdstart = len%4 + 1

    m1 = a[1]
    m2 = a[2]
    m3 = a[3]
    m4 = a[4]
	@inbounds @simd for i in simdstart:4:len
        m1 = ( a[ i ] > m1 ) ? a[ i ] : m1
        m2 = ( a[i+1] > m2 ) ? a[i+1] : m2
        m3 = ( a[i+2] > m3 ) ? a[i+2] : m3
        m4 = ( a[i+3] > m4 ) ? a[i+3] : m4
	end

    return maximum( [ m1, m2, m3, m4 ] )
end


# optimized to use SIMD, performing 4 comparisons at the same time
function maxidx( a::AbstractArray{<:Real,N} ) where {N}

	len = length(a)

	if len < 4
		val, idx = findmax( a ) 
		return idx
	end

	simdstart = len%4 + 1

    maxidx1 = 1
    maxidx2 = 2
    maxidx3 = 3
    maxidx4 = 4
	@inbounds @simd for i in simdstart:4:len
        maxidx1 = ( a[ i ] > a[maxidx1] ) ? i   : maxidx1
        maxidx2 = ( a[i+1] > a[maxidx2] ) ? i+1 : maxidx2
        maxidx3 = ( a[i+2] > a[maxidx3] ) ? i+2 : maxidx3
        maxidx4 = ( a[i+3] > a[maxidx4] ) ? i+3 : maxidx4
	end

    # Finding the max value among the maximum indices computed in parallel
    val, idx = findmax( [  a[ maxidx1 ], a[ maxidx2 ], a[ maxidx3], a[ maxidx4 ] ] );
    maxindex = (idx==1)*maxidx1 + (idx==2)*maxidx2 + (idx==3)*maxidx3 + (idx==4)*maxidx4;

	return maxindex
end

function firstPeak( cmat::Array{<:AbstractFloat,N} ) where {N}

    maxindex = maxidx( cmat )

    # Transforming from linear indexing to cartesian
	h, w, d = size(cmat,1), size(cmat,2), size(cmat,3)

    z = ceil( Int32, maxindex/(h*w) )
    x = ceil( Int32, (maxindex - (z-1)*h*w)/h )
    y = maxindex - (x-1)*h - (z-1)*h*w;

    return (y,x,z), cmat[ maxindex ]
end


function secondPeak( cmat::A{<:AbstractFloat,N}, width::I, idx::NTuple{N,I}) where {N}

	#  First, we store the original values around the maximum peak, and then change them to -1.0.
    ranges = [ max( 1, idx[d] - width ):min( size(cmat,d), idx[d] + width ) for d in 1:N ];
    OGvals = copy( cmat[ ranges... ] );
    cmat[ ranges... ] .= eltype(cmat)(-1.0);

	# Finding the second peak and setting the values back to normal (out of politeness)
    p2, val = firstPeak( cmat );
    cmat[ ranges... ] .= OGvals

    return p2, val
end


setTo0!( A::Array{T,N}... ) where {T<:Number,N} = setTo0!.( A )


function setTo0!( A::Array{T,N} ) where {T<:Number,N}
    zer = T(0)
    @inbounds @simd for i in 1:length(A)
        A[ i ] = zer
    end
end


function SAcoords( IAcoords1, IAcoords2, disp, SM_mp, datasize ) 

	SAcoords1  = max.(    1    , IAcoords1 .+ disp .- SM_mp )
	SAcoords2  = min.( datasize, IAcoords2 .+ disp .+ SM_mp )
	pad_offset = max.( 0, SAcoords1 .- ( IAcoords1 .+ disp .- SM_mp ) )
	y1, x1, z1 = SAcoords1
	y2, x2, z2 = SAcoords2

    return (y1, y2, x1, x2, z1, z2), pad_offset
end


function putWithinSearch!( search::Array{T,N}, data::Array{T,N},
                           so::NTuple{3,I}, coord::NTuple{6,I} 
                         ) where {T<:AbstractFloat,N}

    start  = [ x + 1 for x in so ]; 
    r1, r2 = (coord[1], coord[2])
    c1, c2 = (coord[3], coord[4])
    d1, d2 = (coord[5], coord[6])

	if N == 2
		start[3] = 1; 
		d1, d2 = 1,1;
	end

    @inbounds for z in 0:(d2-d1), c in 0:(c2-c1)
		@simd for r in 0:(r2-r1)
		    search[start[1]+r, start[2]+c, start[3]+z] = data[r1+r, c1+c, d1+z]
		end
    end
end

function putWithinPadded!( padA::Array{Complex{T},N}, A::Array{<:AbstractFloat,N}, mean,
                           offset::NTuple{3,I}, coords::NTuple{6,I}
                         )  where {T<:AbstractFloat,N}

    o1,o2,o3 = offset .+ ( 1, 1, 1 )
    r1, r2 = (coords[1], coords[2])
    c1, c2 = (coords[3], coords[4])
    d1, d2 = (coords[5], coords[6])

    Tzero = T(0)
    Tmean = convert( T, mean )

    for z in 0:(d2-d1), c in 0:(c2-c1)
        @simd for r in 0:(r2-r1)
            padA[ o1+r, o2+c, o3+z ] = complex( convert( T, A[r1+r,c1+c,d1+z] ) - Tmean, Tzero )
        end
    end
end

# construc the integral array of squared values
function integralArraySQ!( padS::AbstractArray{Complex{T},2}, intArr::Array{T,2}
                         ) where {T<:AbstractFloat}

	h, w = size(intArr) .- 1
	@inbounds for c in 2:w+1, r in 2:h+1
		intArr[r,c] = padS[r-1,c-1]^2 + intArr[r-1,c] + intArr[r,c-1] - intArr[r-1,c-1]
	end
end

# optimized sumarea, reusing previous operations
function integralArraySQ!( padS::AbstractArray{Complex{T},3}, intArr::Array{T,3} 
                         ) where {T<:AbstractFloat} 

    h, w, d = size(intArr) .- 1;
    @inbounds for z in 2:d+1, c in 2:w+1
        tmp = 0.0; 
        for r in 2:h+1      
            val2          = real( padS[r-1,c-1,z-1] )^2
            intArr[r,c,z] = val2 + intArr[r,c-1,z] + intArr[r,c,z-1] - intArr[r,c-1,z-1] + tmp;
            tmp          += val2; 
    end end
end

function integralArea( intArr::AbstractArray{<:Real,2}, TL, BR )
	TL   = TL .+ 1;
	BR   = BR .+ 1;
	area = intArr[BR[1],BR[2]] - intArr[BR[1],TL[2]] - intArr[TL[1],BR[2]] + intArr[TL[1],TL[2]]
end

function integralArea( intArr::AbstractArray{<:Real,3}, TL, BR )
	TL = TL .+ 1; 
	BR = BR .+ 1; 
	area  = intArr[BR[1],BR[2],BR[3]] - intArr[TL[1],TL[2],TL[3]]
	area -= intArr[TL[1],BR[2],BR[3]] + intArr[BR[1],TL[2],BR[3]] + intArr[BR[1],BR[2],TL[3]]
	area += intArr[BR[1],TL[2],TL[3]] + intArr[TL[1],BR[2],TL[3]] + intArr[TL[1],TL[2],BR[3]]
    return area
end
