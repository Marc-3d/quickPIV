#=
    Computes and array of magnitudes for 2D or 3D vector fields. Each component of
    the vectors field is provided in a separate array:
        magnitudemap( U, V ) or magnitudemap( U, V, W )
=#

function magnitudemap( arrays::Array{T,N}... ; typ=Float32 ) where {T,N}
    magnitudes = zeros( typ, size( arrays[1] ) )
    magnitudemap!( magnitudes, arrays... )
    return magnitudes
end

function magnitudemap!( magnitudes::Array{U,N}, arrays::Array{T,N}... ) where {T,U,N}
    numcomponents = length( arrays );
    @inbounds for e in 1:length( arrays[1] )
        sum2 = 0.0;
        @simd for c in 1:numcomponents
            sum2 += arrays[c][e]*arrays[c][e]
        end
        magnitudes[e] = convert( U, sqrt( sum2 ) )
    end
end


#=
    Computes a boolean mask, where vectors whose magnitude is > n * std are set to
    true. Each component of the vectors field is provided in a separate array.
=#

function stdmask( n, arrays::Array{T,N}... ) where {T,N}
    mask = zeros( Bool, size( arrays[1] ) )
    stdmask!( n, mask, arrays... );
    return mask
end

function stdmask!( n, mask::Array{Bool,N}, arrays::Array{T,N}... ) where {T,N}
    magnitudes = magnitudemap( arrays... );
    mean = Statistics.mean( magnitudes );
    std  = Statistics.std( magnitudes, mean=mean ) ;
    @simd for e in 1:length( magnitudes )
        @inbounds mask[ e ] = abs( magnitudes[e] - mean ) > n*std
    end
end

#=
    Replacement functions accept a coordinate and a radius argument ( a part from
    the vector field components ). zeroRepalce does not use the coordinate or
    radius, but mean and median require them.
=#
function zeroReplace( coords::NTuple{N,I}, radius::I, arrays::Array{T,N}... ) where {T,N}

    return zeros( T, length(arrays) )
end

function medianReplace( coords::NTuple{N,I}, radius::I, arrays::Array{T,N}... ) where {T,N}

    mincoords = max.( 1, coords .- radius );
    maxcoords = min.( size(arrays[1]), coords .+ radius );
    medians   = zeros( T, length( arrays ) );

    for c in 1:length( arrays )
        medians[ c ] = Statistics.median( arrays[c][ UnitRange.( mincoords, maxcoords )... ] );
    end
    return medians
end

function meanReplace( coords::NTuple{N,I}, radius::I, arrays::Array{T,N}... ) where {T,N}

    mincoords = max.( 1, coords .- radius );
    maxcoords = min.( size(arrays[1]), coords .+ radius );
    means     = zeros( T, length( arrays ) );

    for c in 1:length( arrays )
        means[ c ] = Statistics.mean( arrays[c][ UnitRange.( mincoords, maxcoords )... ] )
    end
    return means
end


"""
    Replaces vectors whose magnitude is > n * std. The replacement function can be
    specified, as well as the radius of the neighbouring vectors that are used to
    compute median and mean vectors.

    Implemented replacement functions are.
        PIV3D.zeroReplace ----> replaces outlier vectors by a vector of zeros
        PIV3D.medianReplace --> replaces outlider vector by the median vector
        PIV3D.meanReplace ----> replaces outlier vectors by the mean vector
"""
function stdFilter( n, arrays::Array{T,N}...; 
				   replace::Function=zeroReplace, radius=0 ) where {T,N}

    filtered = [ copy( arrays[idx] ) for idx in 1:length(arrays) ];

    if n > 0
        mask   = stdmask( n, arrays... );
        coords = CartesianIndices( mask );

        for e in 1:length(mask)
            if mask[e]
                replacement = replace( Tuple(coords[e]), radius, arrays... )
                for c in 1:length(arrays)
                    filtered[c][e] = replacement[c]
                end
            end
        end
    end

    return filtered
end

"""
    Replaces vectors whose magnitude is > n * std. The replacement function can be
    specified, as well as the radius of the neighbouring vectors that are used to
    compute median and mean vectors.

    Implemented replacement functions are.
        PIV3D.zeroReplace ----> replaces outlier vectors by a vector of zeros
        PIV3D.medianReplace --> replaces outlider vector by the median vector
        PIV3D.meanReplace ----> replaces outlier vectors by the mean vector
"""
function magnitudeFilter( threshold, u::Array{T,N}, v::Array{T,N} 
                        ) where {T<:AbstractFloat,N}

    return magnitudeFilter( >, threshold, u, v )
end

function magnitudeFilter( f, threshold, arrays::Array{T,N}...;
                          replace::Function=zeroReplace, radius=0 ) where {T,N}

    magnitudes = magnitudemap( arrays... );
    coords     = CartesianIndices( magnitudes );
    filtered   = [ copy( arrays[idx] ) for idx in 1:length(arrays) ];

    for e in 1:length(magnitudes)
        if f( magnitudes[e], threshold )
            replacement = replace( Tuple(coords[e]), radius, arrays... )
            for c in 1:length(arrays)
                filtered[c][e] = replacement[c]
            end
        end
    end
    return filtered
end


""" Spatial averaging for 2D vector fields"""

function spaceAveraging( avg_radius::I, u::Array{T,2}, v::Array{T,2}; 
                         th=0.1 ) where {T<:AbstractFloat}

    u_avg = zeros( T, size(u) );
    v_avg = zeros( T, size(v) );
    h, w  = size( u )

    for col in 1:w
        cmin = max( 1, col-avg_radius );
        cmax = min( w, col+avg_radius );

        for row in 1:h

            mag = sqrt( u[ row, col ]^2 + v[ row, col ]^2 );
            if mag < th
                continue
            end

            rmin = max( 1, row-avg_radius );
            rmax = min( h, row+avg_radius );
            N = length(cmin:cmax)*length(rmin:rmax);

            mean_u = 0.0;
            mean_v = 0.0;
            @inbounds for x in cmin:cmax
                @simd for y in rmin:rmax
                    mean_u += u[y,x]
                    mean_v += v[y,x]
                end
            end

            u_avg[ row, col ] = mean_u/N;
            v_avg[ row, col ] = mean_v/N;
        end
    end
    return u_avg, v_avg
end

""" Spatial averaging for 3D vector fields"""
function spaceAveraging( avg_radius::I, U::Array{T,3}, V::Array{T,3}, W::Array{T,3} ) where {T<:AbstractFloat}

    u_avg   = zeros( T, size(U) );
    v_avg   = zeros( T, size(V) );
    w_avg   = zeros( T, size(W) );
    h, w, d = size( U )

    for zet in 1:d
        zmin = max( 1, zet-avg_radius );
        zmax = min( d, zet+avg_radius );

        for col in 1:w
            cmin = max( 1, col-avg_radius );
            cmax = min( w, col+avg_radius );

            for row in 1:h

                mag = sqrt( U[row,col,zet]^2 + V[row,col,zet]^2 + W[row,col,zet]^2 )
                if mag == 0
                    continue
                end

                rmin = max( 1, row-avg_radius );
                rmax = min( h, row+avg_radius );
                N = length(rmin:rmax)*length(cmin:cmax)*length(zmin:zmax);

                mean_u = 0.0;
                mean_v = 0.0;
                mean_w = 0.0;

                # Looping over neighbours
                for z in zmin:zmax
                    for x in cmin:cmax
                        for y in rmin:rmax
                            mean_u += U[y,x,z]
                            mean_v += V[y,x,z]
                            mean_w += W[y,x,z]
                        end
                    end
                end

                u_avg[ row, col, zet ] = mean_u/N;
                v_avg[ row, col, zet ] = mean_v/N;
                w_avg[ row, col, zet ] = mean_w/N;
            end
        end
    end

    return u_avg, v_avg, w_avg
end

""" space-time averaging """ # fancier circular mask in notes_16Oct20
function spaceTimeAveraging( avg_rs::I, avg_rt::I, u::Array{T,3}, v::Array{T,3}; th=0.0 )  where {T<:AbstractFloat}

    u_avg = zeros( T, size(u) );
    v_avg = zeros( T, size(v) );

    h, w, time = size( u )

    for t in 1:time
        tmin = max(   1 , t - avg_rt );
        tmax = min( time, t + avg_rt );

        for c in 1:w
            cmin = max( 1, c - avg_rs );
            cmax = min( w, c + avg_rs );

            for r in 1:h

                magnitude = sqrt( u[r,c,t]^2 + v[r,c,t]^2 )
                if  magnitude < th
                    continue
                end

                rmin = max( 1, r - avg_rs );
                rmax = min( h, r + avg_rs );
                N = length(tmin:tmax) + length(rmin:rmax)*length(cmin:cmax) - 1

                # central vector is counted twice, so we account for it by substracting it once
                mean_u = -u[r,c,t];
                mean_v = -v[r,c,t];

                # spatial averaging
                for x in cmin:cmax
                    for y in rmin:rmax
                        mean_u += u[y,x,t]
                        mean_v += v[y,x,t]
                    end
                end

                # temporal averaging
                for tp in tmin:tmax
                    mean_u += u[r,c,tp]
                    mean_v += v[r,c,tp]
                end

                u_avg[r,c,t] = mean_u/N;
                v_avg[r,c,t] = mean_v/N;
            end
        end
    end

    return u_avg, v_avg
end

""" space-time averaging """ # fancier circular mask in notes_16Oct20
function spaceTimeAveraging( avg_rs::I, avg_rt::I, u::Array{T,4}, v::Array{T,4}, w::Array{T,4};
                             th=0.0 )  where {T<:AbstractFloat}

    u_avg = zeros( T, size(u) );
    v_avg = zeros( T, size(v) );
    w_avg = zeros( T, size(w) );

    h, w, d, time = size( u )

    for t in 1:time
        tmin = max(   1 , t - avg_rt );
        tmax = min( time, t + avg_rt );

        for s in 1:d
            smin = max( 1, s - avg_rs );
            smax = min( d, s + avg_rs );

            for c in 1:w
                cmin = max( 1, c - avg_rs );
                cmax = min( w, c + avg_rs );

                for r in 1:h

                    magnitude = sqrt( u[r,c,s,t]^2 + v[r,c,s,t]^2 + w[r,c,s,t]^2 )
                    if  magnitude < th
                        continue
                    end

                    rmin = max( 1, r - avg_rs );
                    rmax = min( h, r + avg_rs );
                    N = length(tmin:tmax) + length(smin:smax)*length(cmin:cmax)*length(rmin:rmax)

                    # central vector is counted twice, so we account for it by substracting it once
                    mean_u = -u[ r, c, s, t ];
                    mean_v = -v[ r, c, s, t ];
                    mean_w = -w[ r, c, s, t ];

                    # spatial averaging, t is constant
                    for z in smin:smax
                        for x in cmin:cmax
                            for y in rmin:rmax
                                mean_u += u[y,x,z,t]
                                mean_v += v[y,x,z,t]
                                mean_w += w[y,x,z,t]
                            end
                        end
                    end

                    # temporal averaging, (r,c,s) are constant
                    for tp in tmin:tmax
                        mean_u += u[r,c,s,tp]
                        mean_v += v[r,c,s,tp]
                        mean_w += w[r,c,s,tp]
                    end

                    u_avg[ r,c,s,t ] = mean_u/N;
                    v_avg[ r,c,s,t ] = mean_v/N;
                    w_avg[ r,c,s,t ] = mean_w/N;
                end
            end
        end
    end

    return u_avg, v_avg, w_avg
end

""" Spatial averaging + similarity thresholding combo """
function similarityAveraging3D( avg_radius::I, U::Array{T,3}, V::Array{T,3}, W::Array{T,3};
                                norm=true, st=0.0 ) where {T<:AbstractFloat}

    u_avg   = zeros( T, size(U) );
    v_avg   = zeros( T, size(V) );
    w_avg   = zeros( T, size(W) );
    h, w, d = size( U )

    for zet in 1:d
        zmin = max( 1, zet-avg_radius );
        zmax = min( d, zet+avg_radius );

        for col in 1:w
            cmin = max( 1, col-avg_radius );
            cmax = min( w, col+avg_radius );

            for row in 1:h

                u1 = U[row,col,zet]
                v1 = V[row,col,zet]
                w1 = W[row,col,zet]
                mag = sqrt( u1*u1 + v1*v1 + w1*w1 )
                if mag == 0
                    continue
                end

                nu1 = U[row,col,zet]/mag
                nv1 = V[row,col,zet]/mag
                nw1 = W[row,col,zet]/mag

                rmin = max( 1, row-avg_radius );
                rmax = min( h, row+avg_radius );

                len  = length(rmin:rmax)*length(cmin:cmax)*length(zmin:zmax);

                # variables to hold the averaged vector component
                mean_u = 0.0;
                mean_v = 0.0;
                mean_w = 0.0;

                # counter of similar neighbouring vectors
                n  = 0;

                for z in zmin:zmax
                    for x in cmin:cmax
                        for y in rmin:rmax

                            mag2 = sqrt( U[y,x,z]^2 + V[y,x,z]^2 + W[y,x,z]^2 )
                            nu2  = U[y,x,z]/mag2
                            nv2  = V[y,x,z]/mag2
                            nw2  = W[y,x,z]/mag2
                            dot  = nu1*nu2 + nv1*nv2 + nw1*nw2

                            if dot > st
                                mean_u += U[y,x,z]
                                mean_v += V[y,x,z]
                                mean_w += W[y,x,z]
                                n += 1;
                            end
                        end
                    end
                end

                sim = (n/len)^2;
                if norm
                    mmag = sqrt( mean_u*mean_u + mean_v*mean_v + mean_w*mean_w );
                else
                    mmag = 1
                end

                u_avg[ row, col, zet ] = mean_u/mmag * sim;
                v_avg[ row, col, zet ] = mean_v/mmag * sim;
                w_avg[ row, col, zet ] = mean_w/mmag * sim;
            end
        end
    end

    return u_avg, v_avg, w_avg
end

""" Spatial averaging """
function similarityMap( avg_radius::I, U::Array{T,3}, V::Array{T,3}, W::Array{T,3};
                        st=0.0 ) where {T<:AbstractFloat}

    similty = zeros( Float32, size(U) );
    h, w, d = size( U )

    for zet in 1:d
        zmin = max( 1, zet-avg_radius );
        zmax = min( d, zet+avg_radius );

        for col in 1:w
            cmin = max( 1, col-avg_radius );
            cmax = min( w, col+avg_radius );

            for row in 1:h

                mag = sqrt( U[row,col,zet]^2 + V[row,col,zet]^2 + W[row,col,zet]^2 )
                if mag == 0
                    continue
                end

                nu1 = U[row,col,zet]/mag
                nv1 = V[row,col,zet]/mag
                nw1 = W[row,col,zet]/mag

                rmin = max( 1, row-avg_radius );
                rmax = min( h, row+avg_radius );

                len  = length(rmin:rmax)*length(cmin:cmax)*length(zmin:zmax);

				len  = avg_radius*avg_radius*avg_radius;

                n = 0;
                for z in zmin:zmax
                    for x in cmin:cmax
                        for y in rmin:rmax
                            mag2 = sqrt( U[y,x,z]^2 + V[y,x,z]^2 + W[y,x,z]^2 )
                            nu2  = U[y,x,z]/mag2
                            nv2  = V[y,x,z]/mag2
                            nw2  = W[y,x,z]/mag2
                            dot  = nu1*nu2 + nv1*nv2 + nw1*nw2
                            if dot > st
                                n += 1;
                            end
                        end
                    end
                end

                similty[row,col,zet] = (n/len);
            end
        end
    end

    return similty
end

function velocityMap( U::Array{T,3}, V::Array{T,3}, W::Array{T,3} ) where {T<:AbstractFloat}

    mags = zeros( Float32, size( U ) )
    for li in 1:length( U )
        mags[li] = sqrt( U[li]*U[li] + V[li]*V[li] + W[li]*W[li] )
    end
    return mags
end

# 2D
function sink( rad::I, pos::II )
    return [ [ (pos[1]-y )/sqrt( (pos[1]-y)^2 + (pos[2]-x)^2 ),
               (pos[2]-x )/sqrt( (pos[1]-y)^2 + (pos[2]-x)^2 ) ] for 
            y in pos[1]-rad:pos[1]+rad, x in pos[2]-rad:pos[2]+rad ]
end

# 3D
function sink( rad::I, pos::III )
    return [ [ (pos[1]-y )/sqrt( (pos[1]-y)^2 + (pos[2]-x)^2 + (pos[3]-z)^2 ),
               (pos[2]-x )/sqrt( (pos[1]-y)^2 + (pos[2]-x)^2 + (pos[3]-z)^2 ),
               (pos[3]-z )/sqrt( (pos[1]-y)^2 + (pos[2]-x)^2 + (pos[3]-z)^2 ) ] for
	       y in pos[1]-rad:pos[1]+rad, x in pos[2]-rad:pos[2]+rad, z in pos[3]-rda:pos[3]+rad ]
end

function crossCorrelateSink!( u, v, w, sink, corr )

    hi, wi, di = size(sink) # interrogation area
    hs, ws, ds = size(u)    # search area

    @inbounds begin
    for z in di:-1:1
        for c in wi:-1:1
            for r in hi:-1:1
                px = sink[ r, c, z ];
                iy = hi-r+1;
                ix = wi-c+1;
                iz = di-z+1;
                corr[ iy:(iy-1+hs),
                      ix:(ix-1+ws),
                      iz:(iz-1+ds) ] .+= px[1] .* u .+ px[2] .* v .+ px[3] .* w
            end
        end
    end
    end # @inbounds
end

function divergenceMap( r::I, U::Array{T,3}, V::Array{T,3}, W::Array{T,3}
                      ) where {T<:AbstractFloat}

    interSink = sink( r, (0,0,0) );

    for idx in 1:length(interSink)
        isnan(interSink[idx][1]) && ( interSink[idx][1] = 0.0 );
        isnan(interSink[idx][2]) && ( interSink[idx][2] = 0.0 );
        isnan(interSink[idx][3]) && ( interSink[idx][3] = 0.0 );
    end

    # normalizing vector fields
    NU = zeros( T, size( U ) );
    NV = zeros( T, size( V ) );
    NW = zeros( T, size( W ) );

    for e in 1:length( U )
        mag   = sqrt( U[e]*U[e] + V[e]*V[e] + W[e]*W[e] )
        NU[e] = ( mag == 0 ) ? 0.0 : U[e]/mag
        NV[e] = ( mag == 0 ) ? 0.0 : V[e]/mag
        NW[e] = ( mag == 0 ) ? 0.0 : W[e]/mag
    end


    cmatrix = zeros( Float32, size(NU) .+ 2*r )

    crossCorrelateSink!( NU, NV, NW, interSink, cmatrix )

    return cmatrix
end


""" PIV Trajectories, work started by Michelle Gottlieb  """

function PIVtrajectories( U::Array{P,4}, V::Array{P,4}, W::Array{P,4},
                          T0, T1, numpoints; facF=1, 
				          subr=( -1: -1, -1: -1, -1: -1 ) ) where {P<:Real}
    return PIVtrajectories( U, V, W, T0, T1, numpoints, 1:size(U,1), 1:size(U,2), 1:size(U,3), subr=subr, facF=facF )
end

function PIVtrajectories( U::Array{P,4}, V::Array{P,4}, W::Array{P,4},
                          T0, T1, numpoints, vrange, lrange, drange; 
                          subr=( -1: -1, -1: -1, -1: -1 ), facF=1 ) where {P<:Real}

    # Vertical, lateral and depth axes.
    VV = collect( vrange );
    LL = collect( lrange );
    DD = collect( drange );

    dims = ( length(VV), length(LL), length(DD) );

    # I am still not sure of what these scaling factors do
    endV = VV[ end ]
    endH = LL[ end ]
    endD = DD[ end ]
    facV = 1 # endV/dims[1]
    facL = 1 # endH/dims[2]
    facD = 1 # endD/dims[3]

    numT = T1 - T0;
    TrajectoriesV = zeros( Float32, numT, numpoints )
    TrajectoriesL = zeros( Float32, numT, numpoints )
    TrajectoriesD = zeros( Float32, numT, numpoints )


    for pidx in 1:numpoints
        # Each iteration in this loop corresponds to a new particle
        # A particle is nothing but a pair of coordinates
        # Random coordinates of the particle inside VV, HH and DD
        c1 = rand( subr[1].start == -1 ? (2:dims[1]-1) : subr[1] )
        c2 = rand( subr[2].start == -1 ? (2:dims[2]-1) : subr[2] )
        c3 = rand( subr[3].start == -1 ? (2:dims[3]-1) : subr[3] )

        # Trajectory starts at the ( VV, HH ) position
        TrajectoriesV[ 1, pidx ] = VV[ c1 ] + 0.0
        TrajectoriesL[ 1, pidx ] = LL[ c2 ] + 0.0
        TrajectoriesD[ 1, pidx ] = DD[ c3 ] + 0.0

        # Translation converted into indices on VV and HH
        Ut = facF * U[ c1, c2, c3, T0 ]
        Vt = facF * V[ c1, c2, c3, T0 ]
        Wt = facF * W[ c1, c2, c3, T0 ]

        # Time-iteration
        for t in 2:numT

            TrajectoriesV[t,pidx] = TrajectoriesV[t-1, pidx] + Ut
            TrajectoriesL[t,pidx] = TrajectoriesL[t-1, pidx] + Vt
            TrajectoriesD[t,pidx] = TrajectoriesD[t-1, pidx] + Wt

            c1 = round( Int32, TrajectoriesV[t,pidx] / facV )
            c2 = round( Int32, TrajectoriesL[t,pidx] / facL )
            c3 = round( Int32, TrajectoriesD[t,pidx] / facD )

            # If new position is out of the coordinate space, stop
            if c1 > dims[1] || c1 < 1 ||
               c2 > dims[2] || c2 < 1 ||
               c3 > dims[3] || c3 < 1
                TrajectoriesV[ t:end, pidx ] .= TrajectoriesV[ t-1, pidx ]
                TrajectoriesL[ t:end, pidx ] .= TrajectoriesL[ t-1, pidx ]
                TrajectoriesD[ t:end, pidx ] .= TrajectoriesD[ t-1, pidx ]
                break
            end

            # Using nextV and nextH to sample the translation from U and V
            Ut = facF * U[ c1, c2, c3, T0+t-1 ]
            Vt = facF * V[ c1, c2, c3, T0+t-1 ]
            Wt = facF * W[ c1, c2, c3, T0+t-1 ]

        end # Time-iteration

    end

    return TrajectoriesV, TrajectoriesL, TrajectoriesD
end
