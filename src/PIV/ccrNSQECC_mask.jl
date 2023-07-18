# include("diffUtils.jl")

"""
	2D IMPLEMENTATION: f is search area, g is interrogation area (translated function)
"""

#  SQERR =  ( IA - SA )² / ( sqrt( sum(IA²) )*sqrt( sum(SA²) ) )
#        = ( IA²      +       SA²     -      2*IA*SA ) / ( sqrt( sum(IA²) )*sqrt( sum(SA²) ) )
#           -----            -----           -------             --------         ---------
#           const     FFTcross(mask,SA^2)  FFTcross(IA,SA)        const      FFTcross(mask,SA^2)

function crossCorrelation!( ::mask_NSQECC, cmat::A{T,2}, shifted::A{T,2}, shifts,
                                           padf::A{C{T},2}, padg::A{C{T},2}, 
                                           padm::A{C{T},2}, padf2::A{C{T},2}, 
                                           plan, iplan
						  ) where {T<:AbstractFloat}

    crossCorrelation!( FFT(), cmat, padf, padg, plan, iplan );  
    crossCorrelation!( FFT(), shifted, padm, padf2, plan, iplan );  

	sumG2 = sum( padg .* padg )
	sqrtsumG2 = sqrt( sumG2 )
	for idx in 1:length(cmat)
		sumF2 = shifted[idx]
		num   = sumG2 + sumF2 - 2*cmat[idx]
		shifted[idx] = num / ( sqrtsumG2 * sqrt( sumF2 ) ); 
	end

    Base.circshift!( cmat, shifted, shifts[1:2] ) # circshift!( dest, src, shifts )
end


function PIV_2D( ::mask_NSQECC, img1::A{<:Real,2}, img2::A{<:Real,2}, mask::A{<:Bool,2},
                                IA::III, SM::III, overlap::III, mpass::I, width::I, 
                                peak::S, sigNoise::S, filtFun::Function, threshold::F ; 
                                corrType=Float32, vfType=Float32 ) 
    
	# Calculating size of the vector field
    step   = IA .- overlap; 
	imsize = ( size(img1,1), size(img1,2), size(img1,3) );
    VFsize = length.( StepRange.( IA, step, imsize ) ); 

    U  = zeros( vfType, VFsize[1:2] ); 
    V  = zeros( vfType, VFsize[1:2] ); 
    SN = zeros( vfType, VFsize[1:2] );

	ignoreSN = ( sigNoise == "" ); 
	SM = ( SM[1], SM[2], 0 ); 
    
    for mp in mpass:-1:1

        last_mp = ( mp == 1 ); 
                
        # Scaling IA, SM and step (overlap) to multi-pass iteration.
        IA_mp    =  IA  .* mp; 
        SM_mp    =  SM  .* mp; 
        SA_mp    = IA_mp .+ 2 .* SM_mp
        step_mp  = step .* mp;
		IAranges = StepRange.( IA_mp, step_mp, imsize )  

        # Initialize cross-correlation variables
        csize   = 2 .* ( IA_mp .+ SM_mp  );
        shifts  = div.( csize, 2 ) .+ SM_mp .- 1; 
        cmatrix = zeros( corrType, csize[1:2] );
        shifted = zeros( corrType, csize[1:2] ); 
        pads    = zeros( Complex{corrType}, csize[1:2] ); 
        padi    = zeros( Complex{corrType}, csize[1:2] );
		padm    = zeros( Complex{corrType}, csize[1:2] );
		pads2   = zeros( Complex{corrType}, csize[1:2] );
        plan    =  FFTW.plan_fft!( pads ); 
        iplan   = FFTW.plan_ifft!( pads );

		# Progress bar
		n = length( IAranges[2] ) * length( IAranges[1] )
		p = ProgressMeter.Progress(n, "Computing PIV...")

		# Nested loop establishing IV and UVW coordinates
		vfx = [ 1-mp, 0 ]
		for x2 in IAranges[2];     x1 = x2-IA_mp[2]+1; vfx .+= mp; vfy = [ 1-mp, 0 ];
			for y2 in IAranges[1]; y1 = y2-IA_mp[1]+1; vfy .+= mp;
				
				ProgressMeter.next!( p )

				# 0-. Filtering undesirable IAs, ex background
				interr = view( img1, y1:y2, x1:x2 )

				if filtFun( interr ) < threshold 
					SN[ vfy[1]:vfy[2], vfx[1]:vfx[2] ] .= NaN32
					continue; 
				end

			    # 1.1-. Previously computed displacements shift the search area during multipass
			    u = round(Int32, U[ vfy[1], vfx[1] ]);
			    v = round(Int32, V[ vfy[1], vfx[1] ]);

				scoords, offset = SAcoords( (y1,x1,1), (y2,x2,1), (u,v,0), SM_mp, imsize );
				    
				# 1.2-. Copying IV and SV into the padded arrays for FFT
				setTo0!( pads )
				setTo0!( padi )
				setTo0!( padm ); 
				setTo0!( pads2 ); 
			    putWithinPadded!( pads, img2, 0.0, offset, scoords );
			    putWithinPadded!( pads2, img2 .^ 2, 0.0, offset, scoords );
			    putWithinPadded!( padi, img1, 0.0, (0,0,0), (y1,y2,x1,x2,1,1) );
				padi[ y1:y2, x1:x2 ] .*= mask[ y1:y2, x1:x2 ];
			    putWithinPadded!( padm, mask, 0.0, (0,0,0), (y1,y2,x1,x2,1,1) );

				# 2-. Cross-correlation
				crossCorrelation!( mask_NSQECC(), cmatrix, shifted, shifts, pads, padi, padm, pads2, plan, iplan, SA_mp, IA_mp )

				# 3.1-. Calculation of displacement
		        ( r, c ) = approxTranslation( cmatrix, peak, Val(last_mp) )

				# 3.2-. Computing signal to noise ration, if needed
				if ( !ignoreSN && last_mp )
		        	SN[ vfy[1]:vfy[2], vfx[1]:vfx[2] ] = SNRatio( cmatrix, sigNoise, width )
				end

		        # 4-. Updating U, V matrices 
		        U[ vfy[1]:vfy[2], vfx[1]:vfx[2] ] .= u - r;
		        V[ vfy[1]:vfy[2], vfx[1]:vfx[2] ] .= v - c;
		end	end
	end # mulitpass
    
    return U, V, SN
end


"""
	3D IMPLEMENTATION: f is search volume, g is interrogation volume (translated function)
"""

function crossCorrelation!( ::NSQECC, cmat::A{T,3}, shifted::A{T,3}, shifts,
                                      padF::A{C{T},3}, padG::A{C{T},3}, plan, iplan,
                                      intArr2::A{T,3}, sumG2, sizeF, sizeG
						  ) where {T<:AbstractFloat}

	# Computing IV*SV
    crossCorrelation!( FFT(), cmat, padF, padG, plan, iplan );  
    Base.circshift!( shifted, cmat, shifts )
	
	# Computing sqrt(sum(IV²))
	denG = sqrt( sumG2 )

	# Computing sum(SV²) and sqrt(sum(SV²)). 
	fh, fw, fd = sizeF
	gh, gw, gd = sizeG .- 1

    for z in 0:size(cmat,3)-1;         z1 = max( 1, fd - z ); z2 = min( fd, z1 + gd )
		for c in 0:size(cmat,2)-1;     c1 = max( 1, fw - c ); c2 = min( fw, c1 + gw )
			for r in 0:size(cmat,1)-1; r1 = max( 1, fh - r ); r2 = min( fh, r1 + gh )

				sumF2 = abs( integralArea( intArr2, (r1-1,c1-1,z1-1), (r2,c2,z2) ) )
				num   = sumF2 + sumG2 - 2*shifted[r+1,c+1,z+1]

				cmat[r+1,c+1,z+1] = 1/( 1 + num/( denG * sqrt(sumF2) ) )
	end	end	end
end

function PIV_3D( ::NSQECC, vol1::A{<:Real,3}, vol2::A{<:Real,3},
                         IA::III, SM::III, overlap::III, mpass::I, width::I, 
						 peak::S, sigNoise::S, filtFun::Function, threshold::F ;
                         corrType=Float32, vfType=Float32 ) 
    
	# Calculating size of the vector field
    step   = IA .- overlap; 
    VFsize = length.( StepRange.( IA, step, size(vol1) ) ); 

    U  = zeros( vfType, VFsize ); 
    V  = zeros( vfType, VFsize ); 
    W  = zeros( vfType, VFsize );
    SN = zeros( vfType, VFsize );
	df = zeros( vfType, VFsize ); 

	ignoreSN = ( sigNoise == "" ); 
    
    for mp in mpass:-1:1

        last_mp = ( mp == 1 ); 
                
        # Scaling IA, SM and step (overlap) to multi-pass iteration.
        IA_mp   =  IA  .* mp; 
        SM_mp   =  SM  .* mp; 
        step_mp = step .* mp;
		IAranges = StepRange.( IA_mp, step_mp, size(vol1) )

        # Initialize cross-correlation variables
        csize   = 2 .* ( IA_mp .+ SM_mp  );
        cmatrix = zeros( corrType, csize );
        shifted = zeros( corrType, csize ); 
        shifts  = div.(csize,2) .+ SM_mp .- 1; 
        pads    = zeros( Complex{corrType}, csize ); 
        padi    = zeros( Complex{corrType}, csize ); 
        plan    =  FFTW.plan_fft!( pads ); 
        iplan   = FFTW.plan_ifft!( pads );

		# Initializing integral array for search area
        ssize   = IA_mp .+ 2 .* SM_mp
        intArr2 = zeros( corrType, ssize .+ 1 ); 

		n = length( IAranges[3] ) * length( IAranges[2] ) * length( IAranges[1] )
		p = ProgressMeter.Progress(n, "Computing PIV...")

		# Nested loop establishing IV and UVW coordinates
		vfz = [ 1-mp, 0 ]
		for z2 in IAranges[3];         z1 = z2-IA_mp[3]+1; vfz .+= mp; vfx = [ 1-mp, 0 ];
			for x2 in IAranges[2];     x1 = x2-IA_mp[2]+1; vfx .+= mp; vfy = [ 1-mp, 0 ];
				for y2 in IAranges[1]; y1 = y2-IA_mp[1]+1; vfy .+= mp;

					ProgressMeter.next!( p )
				
					interr = view( vol1, y1:y2, x1:x2, z1:z2 )
					
					if filtFun( interr ) < threshold # Filtering of undesirable IAs, ex background
						SN[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= -1.0
						continue; 
					end

				    # 1-. Previously computed displacements shift the search area during multipass
				    u = round(Int32, U[ vfy[1], vfx[1], vfz[1] ]);
				    v = round(Int32, V[ vfy[1], vfx[1], vfz[1] ]);
				    w = round(Int32, W[ vfy[1], vfx[1], vfz[1] ]);

				    # 1.2-. Search volume coordinates after shifting.
					scoords, so = SAcoords( (y1,x1,z1), (y2,x2,z2), (u,v,w), SM_mp, size(vol1) );

					#println( (y1,y2,x1,x2,z1,z2), scoords, so, vfy, vfx, vfz, (u,v,w) )
				    
					# 1.3-. Copying IV and SV into the padded arrays for FFT
					setTo0!( pads )
					setTo0!( padi )
				    putWithinPadded!( pads, vol2, 0.0, so, scoords );
				    putWithinPadded!( padi, vol1, 0.0, (0,0,0), (y1,y2,x1,x2,z1,z2) );

					# 1.4-. Constructing the integral array of SV
					setTo0!( intArr2 )
					integralArraySQ!( pads, intArr2 )

					# 1.5 -. Power of IA and SA
					sumG2 = 0.0
					@inbounds for t3 in z1:z2, t2 in x1:x2
						@simd for t1 in y1:y2
							sumG2 += vol1[t1,t2,t3]*vol1[t1,t2,t3]
					end	end

					# 2-. Cross-correlation
					crossCorrelation!( NSQECC(), cmatrix, shifted, shifts, pads, padi, plan, iplan, 
                                               intArr2, sumG2, ssize, IA_mp )

				    # 3-. Extracting the displacement
				    ( r, c, z ) = approxTranslation( cmatrix, peak, Val(last_mp) )

					# 3.2 -. Computing signal to noise ration, if needed
		 			if ( !ignoreSN && last_mp )
						snres = SNRatio( cmatrix, sigNoise, width )
				    	SN[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= snres
					end

				    # 4-. Updating U, V, W matrices 

				    U[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= u - r;
				    V[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= v - c;
				    W[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= w - z;
				   df[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= maxval( cmatrix );
		end	end	end
	end # multipass
    
    return U, V, W, SN, df
end
