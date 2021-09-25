"""
	2D IMPLEMENTATION
"""

# crosscorr( f, g ) = iFFT( conj( FFT( f ) ) .* FFT( g ) )

function crossCorrelation!( ::FFT, cmat::A{T,2}, padf::A{C{T},2}, padg::A{C{T},2}, plan, iplan
                          ) where {T<:AbstractFloat} 

    plan * padf;
	plan * padg;
    @inbounds @simd for e in 1:length(padf)
        padf[e] = conj( padf[e] ) * padg[e]; 
    end
    
	iplan * padf; 
	@inbounds @simd for e in 1:length(padf) 
		cmat[e] = real( padf[e] ) 
	end
end

function PIV_2D( ::FFT, img1::A{<:Real,2}, img2::A{<:Real,2},
                        IA::III, SM::III, overlap::III, mpass::I, width::I,
						peak::S, sigNoise::S, filtFun::Function, threshold::F;
 					    corrType=Float32, vfType=Float32 ) 
        
    step   = IA  .-  overlap; 
	imsize = ( size(img1,1), size(img1,2), size(img1,3) )
    vfSize = length.( StepRange.( IA, step, imsize ) ); 

    U  = zeros( vfType, vfSize[1:2] ); 
    V  = zeros( vfType, vfSize[1:2] ); 
    SN = zeros( vfType, vfSize[1:2] );
	
	ignoreSN = sigNoise == "";
	SM = ( SM[1], SM[2], 0 );  

    for mp in mpass:-1:1

        last_mp = ( mp == 1 ); 
                  
        # Scaling IA, SM and step (overlap) to multi-pass.
        IA_mp    =  IA  .* mp; 
        SM_mp    =  SM  .* mp; 
        step_mp  = step .* mp; 
		IAranges = StepRange.( IA_mp, step_mp, imsize ) 
                       
        # Initializing variables for FFT cross-Correlation. 
        csize   = 2 .* ( IA_mp .+ SM_mp )
        shifts  = div.( csize, 2 ) .+ SM_mp .- 1;
        cmatrix = zeros( corrType, csize[1:2] );
        shifted = zeros( corrType, csize[1:2] ); 
        pads    = zeros( Complex{corrType}, csize[1:2] );
        padi    = zeros( Complex{corrType}, csize[1:2] ); 
        plan    = FFTW.plan_fft!(  pads ); 
        iplan   = FFTW.plan_ifft!( pads ); 

		# Nested loop establishing IA and UV coordinates
		vfx = [ 1-mp, 0 ]
		for x2 in IAranges[2];     x1 = x2-IA_mp[2]+1; vfx .+= mp; vfy = [ 1-mp, 0 ];
			for y2 in IAranges[1]; y1 = y2-IA_mp[1]+1; vfy .+= mp;

				# 0-. Filtering undesirable IAs, ex background
				interr = view( img1, y1:y2, x1:x2 )

				if filtFun( interr ) < threshold 
					SN[ vfy[1]:vfy[2], vfx[1]:vfx[2] ] .= -1.0
					continue; 
				end

			    # 1.1-. Previously computed displacements shift the search area during multipass
			    u = round(Int32, U[ vfy[1], vfx[1] ]);
			    v = round(Int32, V[ vfy[1], vfx[1] ]);

				# 1.2-. Search volume coordinates after shifting.
				scoords, offset = SAcoords( (y1,x1,1), (y2,x2,1), (u,v,0), SM_mp, imsize );

				# 1.3-. Copying IA and SA into the padded arrays for FFT
				setTo0!( pads, padi )
			    putWithinPadded!( pads, img2, 0.0, offset, scoords );
			    putWithinPadded!( padi, img1, 0.0, (0,0,0), (y1,y2,x1,x2,1,1) );

           		# 2-. cross-correlation
		        crossCorrelation!( FFT(), cmatrix, pads, padi, plan, iplan )
		        Base.circshift!( shifted, cmatrix, shifts[1:2] )

			    # 3.1-. Extracting the displacement
			    ( r, c ) = approxTranslation( shifted, peak, Val(last_mp) )

				# 3.2 -. Computing signal to noise ration, if needed
	 			if ( !ignoreSN && last_mp )
			    	SN[ vfy[1]:vfy[2], vfx[1]:vfx[2] ] .= SNRatio( shifted, sigNoise, width )
				end
				
			    # 4-. Updating U, V matrices 
			    U[ vfy[1]:vfy[2], vfx[1]:vfx[2] ] .= u - r;
			    V[ vfy[1]:vfy[2], vfx[1]:vfx[2] ] .= v - c;
		end end
    end # Multi-pass loop

    return U, V, SN
end

"""
    3D IMPLEMENTATION
"""

# crosscorr( f, g ) = iFFT( conj( FFT( f ) ) .* FFT( g ) )
function crossCorrelation!( ::FFT, cmat::A{T,3}, padf::A{C{T},3}, padg::A{C{T},3}, plan, iplan 
                          ) where {T<:AbstractFloat} 
    plan * padf;
    plan * padg;
    @inbounds @simd for e in 1:length(padf)
        padf[e] = conj( padf[e] ) * padg[e]; 
    end

    iplan * padf; 
    @inbounds @simd for e in 1:length(padf) 
        cmat[e] = real( padf[e] ) 
    end
end

function PIV_3D( ::FFT, vol1::A{T,3}, vol2::A{T,3},
                        IA::III, SM::III, overlap::III, mpass::I, width::I,
                        peak::S, sigNoise::S, filtFun::Function, threshold::F ; 
                        corrType=Float32, vfType=Float32
               ) where {T<:Real} 
    
	# Calculating vector field size.
    step   = IA  .-  overlap; 
    vfSize = length.( StepRange.( IA, step, size(vol1) ) ); 

    U  = zeros( vfType, vfSize ); 
    V  = zeros( vfType, vfSize ); 
    W  = zeros( vfType, vfSize );
    SN = zeros( vfType, vfSize );

	ignoreSN = sigNoise == ""; 

    for mp in mpass:-1:1

        last_mp = ( mp == 1 ); 

        # Scaling IA, SM and step (overlap) to multi-pass iteration.
        IA_mp    =  IA  .* mp; 
        SM_mp    =  SM  .* mp; 
        step_mp  = step .* mp; 
		IAranges = StepRange.( IA_mp, step_mp, size(vol1) ) 
            
        # Initialize cross-correlation variables once per multi-pass iteration
        csize   = 2 .* ( IA_mp .+ SM_mp )
        shifts  = div.( csize, 2 ) .+ SM_mp .- 1; 
        cmatrix = zeros( corrType, csize );
        shifted = zeros( corrType, csize ); 
        padi    = zeros( Complex{corrType}, csize ); 
        pads    = zeros( Complex{corrType}, csize );
        plan    = FFTW.plan_fft!(  pads ); 
        iplan   = FFTW.plan_ifft!( pads );

		# Nested loop establishing IV and UVW coordinates
		vfz = [ 1-mp, 0 ]
		for z2 in IAranges[3];         z1 = z2-IA_mp[3]+1; vfz .+= mp; vfx = [ 1-mp, 0 ];
			for x2 in IAranges[2];     x1 = x2-IA_mp[2]+1; vfx .+= mp; vfy = [ 1-mp, 0 ];
				for y2 in IAranges[1]; y1 = y2-IA_mp[1]+1; vfy .+= mp;

					# 0-. Filtering undesirable IAs, ex background
					interr = view( vol1, y1:y2, x1:x2, z1:z2 )
					
					if filtFun( interr ) < threshold
						SN[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= NaN32
						continue; 
					end

				    # 1-. Previously computed displacements shift the search area during multipass
				    u = round(Int32, U[ vfy[1], vfx[1], vfz[1] ]);
				    v = round(Int32, V[ vfy[1], vfx[1], vfz[1] ]);
				    w = round(Int32, W[ vfy[1], vfx[1], vfz[1] ]);

				    # 1.2-. Search volume coordinates after shifting.
					scoords, off = SAcoords( (y1,x1,z1), (y2,x2,z2), (u,v,w), SM_mp, size(vol1) );
				    
					# 1.3-. Copying IV and SV into the padded arrays for FFT
					setTo0!( pads )
					setTo0!( padi )
				    putWithinPadded!( pads, vol2, 0.0, off, scoords );
				    putWithinPadded!( padi, vol1, 0.0, (0,0,0), (y1,y2,x1,x2,z1,z2) );

					# 2-. Cross-correlation
					crossCorrelation!( FFT(), cmatrix, pads, padi, plan, iplan )
		            Base.circshift!( shifted, cmatrix, shifts );

				    # 3-. Extracting the displacement
				    ( r, c, z ) = approxTranslation( shifted, peak, Val(last_mp) )

					# 3.2 -. Computing signal to noise ration, if needed
		 			if ( !ignoreSN && last_mp )
						snres = SNRatio( shifted, sigNoise, width )
				    	SN[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= snres
					end

				    # 4-. Updating U, V, W matrices 
				    U[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= u - r;
				    V[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= v - c;
				    W[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= w - z;
		end	end	end             
    end # Multi-pass loop
    
    return U, V, W, SN
end



""" FFT out of place CROSS CORRELATION """

function crossCorrelation( ::FFT, g::A{<:Real,N}, f::A{<:Real,N}; corrType=Float32, shift=false ) where {N} 
    
    csize = size( f ) .+ size( g ); 

    corr  = zeros( corrType, csize ); 
    padf  = zeros( Complex{corrType}, csize ); 
    padg  = zeros( Complex{corrType}, csize );
    plan  =  plan_fft!( padf );
    iplan = plan_ifft!( padf ); 

    putWithinPadded!( padf, f, 0.0, (0,0,0), (1,size(f,1),1,size(f,2),1,size(f,3)) ); 
    putWithinPadded!( padg, g, 0.0, (0,0,0), (1,size(g,1),1,size(g,2),1,size(f,3)) ); 

	crossCorrelation!( FFT(), corr, padf, padg, plan, iplan )
   
    if shift
        shifted = copy( corr ); 
        shifts  = div.( csize, 2 ) .+ div.( size(f) .- size(g), 2 ) .- 1; 
        Base.circshift!( corr, shifted, shifts );
    end
    
   return corr
end
