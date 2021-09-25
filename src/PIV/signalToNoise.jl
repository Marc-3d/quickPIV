"""
   SIGNAL TO NOISE MEASUREMENT
"""
function SNRatio( corr::Array{T,N}, method::S, width::I ) where {T<:Real,N}
        
    max2       = 0.0; 
    h, w, d    = size(corr,1), size(corr,2), size(corr,3)
    idx1, max1 = firstPeak( corr )
    
    if     method == "peak2peak"
		#idx2, max2 = secondPeak( corr, width, idx1[1:N] )
        idx2, max2 = value2ndPeak( corr, max1, width )
    elseif method == "peak2mean"
        max2 = Statistics.mean( corr )
    else
        throw(ErrorException("Error: use 'peak2peak' or 'peak2mean'.\n"))
    end

    # if max2 == 0, sig2noise = inf. Check this at return

	# When using ZNCC, max2 might be negative, then max1/max2 is negative value
	# and we can't take the logarithm. Therefore, shift max1 and max2 up, until 
	# max2 = 1 and max1 = max1 - max2 + 1. Then we can divide max1/max2 without problems
	if max2 < 0 && max1 > 0 
		max1 = max1 - max2 + 1
		max2 = 1 # max2 - max2 + 1
	end

    sig2noise = max2 !== 0.0 ? max1 / max2 : 0.0

    return sig2noise; 
end

function value2ndPeak( cmat::Array{T,N}, max1, rad ) where {T<:AbstractFloat,N}
    
	h, w, d = size(cmat,1), size(cmat,2), size(cmat,3); 
	zrad    = ( N == 2 ) ? 0 : rad; 
    max2    = cmat[ 1, 1, 1 ]
    maxid   = ( 1, 1, 1 )
    
    @inbounds for dep in 1+zrad:d-zrad, col in 1+rad:w-rad, row in 1+rad:h-rad
        
        ( cmat[row,col,dep] == max1 || cmat[row,col,dep] <= cmat[maxid...] ) && ( continue; );
            
        ismax = true; 
        for z in dep-zrad:dep+zrad, c in col-rad:col+rad, r in row-rad:row+rad
            ( cmat[row,col,dep] < cmat[r,c,z]  ) && ( ismax = false; break; )
        end
        
        max2  = ( ismax ) ? cmat[row,col,dep] : max2; 
		#maxid = ( ismax ) ? ( row, col, dep ) : maxid;  
    end
    
    return maxid, max2
end
