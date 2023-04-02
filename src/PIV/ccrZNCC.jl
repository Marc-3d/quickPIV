# include("znccUtils.jl")

"""
	2D IMPLEMENTATION: f is search volume, g is interrogation volume (translated function) 
"""

#  Lewis, J.P.: Fast normalized cross-correlation. Ind. Light Magic10(2001)
#
#  ZNCC = (IA - meanIA )(SA - meanSA ) / ( sqrt( sum(IA-meanIA)^2 )*sqrt( sum(SA-meanSA)^2 ) )
#          --------------------------            ----------------         ----------------
#        FFTcross( IA-meanIA, SA-meanSA)    constant: std(IA)*sqrt(N)       integralArray

function crossCorrelation!( ::ZNCC, cmat::A{T,2}, shifted::A{T,2}, shifts,
                                    padf::A{C{T},2}, padg::A{C{T},2}, plan, iplan, 
                                    intArr2::A{T,2}, stdG, sizeF, sizeG
                           ) where {T<:AbstractFloat}
    
    # Computing numerator 
    crossCorrelation!( FFT(), cmat, padf, padg, plan, iplan );
    Base.circshift!( shifted, cmat, shifts[1:2] ) 
    
	# Constant denominator
    N    = prod(sizeG); 
	denG = stdG*sqrt(N);

	# Variable denominator for each translation
	fh, fw = sizeF
	gh, gw = sizeG .- 1

	for c in 0:size(cmat,2)-1; 		c1 = max(1, fw - c ); c2 = min( fw, c1 + gw ) 
		for r in 0:size(cmat,1)-1;  r1 = max(1, fh - r ); r2 = min( fh, r1 + gh )

			sumF2 = integralArea( intArr2, (r1-1,c1-1), (r2,c2) ) 
			sumF2 = ( sumF2 <= 0.0 ) ? Inf : sumF2; 
		    cmat[ r+1, c+1 ] = shifted[ r+1, c+1 ]/( sqrt( sumF2 )*denG );
    end	end
end


function PIV_2D( ::ZNCC, img1::A{<:Real,2}, img2::A{<:Real,2},
                         IA::III, SM::III, overlap::III, mpass::I, width::I,
                         peak::S, sigNoise::S, filtFun::Function, threshold::F; 
					     corrType=Float32, vfType=Float32 )
    
    step   = IA  .-  overlap; 
	imsize = ( size(img1,1), size(img1,2), size(img1,3) )
    vfSize = length.( StepRange.( IA, step, imsize ) );
 
    U  = zeros( vfType, vfSize[1:2] ); 
    V  = zeros( vfType, vfSize[1:2] ); 
    SN = zeros( vfType, vfSize[1:2] );

	ignoreSN = ( sigNoise == "" );
	SM = ( SM[1], SM[2], 0 );  
    
    for mp in mpass:-1:1

        last_mp = ( mp == 1 ); 
                
        # Scaling IA, SM and step (overlap) to multi-pass iteration.
        IA_mp    =  IA  .* mp; 
        SM_mp    =  SM  .* mp; 
        step_mp  = step .* mp;
		IAranges = StepRange.( IA_mp, step_mp, imsize )  

        # Initializing variables for FFT cross-Correlation. 
        csize   = 2 .* ( IA_mp .+ SM_mp  );
		shifts  = div.(csize , 2) .+ SM_mp .- 1; 
        cmatrix = zeros( corrType, csize[1:2] );
        num     = zeros( corrType, csize[1:2] ); 
        pads    = zeros( Complex{corrType}, csize[1:2] ); 
        padi    = zeros( Complex{corrType}, csize[1:2] ); 
        plan    =  FFTW.plan_fft!( pads ); 
        iplan   = FFTW.plan_ifft!( pads );

		# Initializing integral arrays for computing ZNCC
        ssize   = IA_mp .+ 2 .* SM_mp
        intArr2 = zeros( corrType, ssize[1:2] .+ 1 ); 

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

				scoords, offset = SAcoords( (y1,x1,1), (y2,x2,1), (u,v,0), SM_mp, imsize );

				# 1.3-. Computing the mean and std of IA and mean and mean^2 of SA
				meanI, stdI = meanStd( interr )

				meanS = 0.0
				@inbounds for x in x1:x2
					@simd for y in y1:y2
						meanS += img2[y,x]
				end	end
				meanS = meanS/length(img2)

				# 1.4-. Copying IA - meanI and SA - meanS into the padded arrays for FFT
				setTo0!( pads, padi )
			    putWithinPadded!( pads, img2, meanS, offset, scoords );
			    putWithinPadded!( padi, img1, meanI, (0,0,0), (y1,y2,x1,x2,1,1) );

				# 1.5-. Computing integral array and integral array ^2
				setTo0!( intArr2 )
				integralArraySQ!( pads, intArr2 )

		        # 2-. Cross-Correlation
		        crossCorrelation!( ZNCC(), cmatrix, num, shifts, pads, padi, plan, iplan, 
		                                   intArr2, stdI, ssize, IA_mp );

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
    end # Multi-pass loop
    
    return U, V, SN
end
    
"""
	3D IMPLEMENTATION: f is search volume, g is interrogation volume (translated function)
"""

function crossCorrelation!( ::ZNCC, cmat::A{T,3}, shifted::A{T,3},  shifts,
                                    padf::A{C{T},3}, padg::A{C{T},3}, plan, iplan,
                                    intArr2::A{T,3}, stdG, sizeF, sizeG
						  ) where {T<:AbstractFloat}
    
    # Computing numerator 
    crossCorrelation!( FFT(), cmat, padf, padg, plan, iplan );  
    Base.circshift!( shifted, cmat, shifts )

	# Constant denominator
	N = prod(sizeG); 
	denG = stdG*sqrt(N);

	# Variable denominator for each translation
	fh, fw, fd = sizeF
	gh, gw, gd = sizeG .- 1
	off        = sizeF .- sizeG .+ 1 # = 2*searchMargin + 1

    for z in 0:size(cmat,3)-1;         #z1 = max( 1, fd - z ); z2 = min( fd, z1 + gd )
                                       z3 = max( 1 ,min(off[3],fd-z)); z4 = z3 + gd;  
		for c in 0:size(cmat,2)-1;     #c1 = max( 1, fw - c ); c2 = min( fw, c1 + gw )
                                       c3 = max( 1 ,min(off[2],fw-c)); c4 = c3 + gw;  
			for r in 0:size(cmat,1)-1; #r1 = max( 1, fh - r ); r2 = min( fh, r1 + gh )
									   r3 = max( 1 ,min(off[1],fh-r)); r4 = r3 + gh;

				sumF2 = integralArea( intArr2, (r3-1,c3-1,z3-1), (r4,c4,z4) )
				sumF2 = ( sumF2 <= 0.0 ) ? Inf : sumF2; 
                cmat[ r+1, c+1, z+1 ] = shifted[ r+1, c+1, z+1 ]/( sqrt( sumF2 )*denG );
	end	end end
end

function PIV_3D( ::ZNCC, vol1::A{T,3}, vol2::A{T,3},
                         IA::III, SM::III, overlap::III, mpass::I, width::I,
                         peak::S, sigNoise::S, filtFun::Function, threshold::F ; 
                         corrType=Float32, vfType=Float32
			   ) where {T<:Real} 
    
    step   = IA  .-  overlap; 
    vfSize = length.( StepRange.( IA, step, size(vol1) ) ); 

    U  = zeros( vfType, vfSize ); 
    V  = zeros( vfType, vfSize ); 
    W  = zeros( vfType, vfSize );
    SN = zeros( vfType, vfSize );
	df = zeros( vfType, vfSize ); 

	ignoreSN = sigNoise == ""; 

	p = ProgressMeter.Progress(n, dt=1.0)
    
    @inbounds for mp in mpass:-1:1
        last_mp = ( mp == 1 ); 
                
        # Scaling IA, SM and step (overlap) to multi-pass iteration.
        IA_mp    =  IA  .* mp; 
        SM_mp    =  SM  .* mp; 
        step_mp  = step .* mp; 
		IAranges = StepRange.( IA_mp, step_mp, size(vol1) ) 
   
		# Initializing integral arrays for computing ZNCC
        ssize   = IA_mp .+ 2 .* SM_mp
        intArr2 = zeros( corrType, ssize .+ 1 ); 

        # Initialize cross-correlation variables once per multi-pass iteration
        csize   = 2 .* ( IA_mp .+ SM_mp  );
        cmatrix = zeros( corrType, csize );
        shifted = zeros( corrType, csize ); 
        shifts  = div.(csize,2) .+ SM_mp .- 1; 
        pads    = zeros( Complex{corrType}, csize ); 
        padi    = zeros( Complex{corrType}, csize ); 
        plan    =  FFTW.plan_fft!( pads ); 
        iplan   = FFTW.plan_ifft!( pads );

		# Nested loop establishing IV and UVW coordinates
		vfz = [ 1-mp, 0 ]
		for z2 in IAranges[3];         z1 = z2-IA_mp[3]+1; vfz .+= mp; vfx = [ 1-mp, 0 ];
			for x2 in IAranges[2];     x1 = x2-IA_mp[2]+1; vfx .+= mp; vfy = [ 1-mp, 0 ];
				for y2 in IAranges[1]; y1 = y2-IA_mp[1]+1; vfy .+= mp;

					ProgressMeter.next!( p )

					# 0-. Filtering undesirable IAs, ex background
					interr = view( vol1, y1:y2, x1:x2, z1:z2 )
					
					if filtFun( interr ) < threshold
						SN[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= -1.0
						continue; 
					end

				    # 1-. Previously computed displacements shift the search area during multipass
				    u = round(Int32, U[ vfy[1], vfx[1], vfz[1] ]);
				    v = round(Int32, V[ vfy[1], vfx[1], vfz[1] ]);
				    w = round(Int32, W[ vfy[1], vfx[1], vfz[1] ]);

					# 1.2-. Search volume coordinates after shifting.
					scoords, off = SAcoords( (y1,x1,z1), (y2,x2,z2), (u,v,w), SM_mp, size(vol1) );
				    
					# 1.3-. Computing the mean and std of IA and mean and mean^2 of SA
					meanI, stdI = meanStd( interr )

					meanS = 0.0
					@inbounds for z in z1:z2, x in x1:x2
						@simd for y in y1:y2
							meanS += vol2[y,x,z]
					end	end
					meanS = meanS/length(vol2)

					# 1.4-. Copying IA - meanI and SA - meanS into the padded arrays for FFT
					setTo0!( pads, padi )
					putWithinPadded!( pads, vol2, meanS, off, scoords );
					putWithinPadded!( padi, vol1, meanI, (0,0,0), (y1,y2,x1,x2,z1,z2) );

					# 1.5-. Computing integral array and integral array ^2
					setTo0!( intArr2 )
					integralArraySQ!( pads, intArr2 )

				    # 2-. Cross-Correlation
				    crossCorrelation!( ZNCC(), cmatrix, shifted, shifts, pads, padi, plan, iplan, 
				                       intArr2, stdI, ssize, IA_mp );

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
				   df[ vfy[1]:vfy[2], vfx[1]:vfx[2], vfz[1]:vfz[2] ] .= maxval(cmatrix)
		end	end	end
    end # Multi-pass loop
    
    return U, V, W, SN, df
end


""" ZNCC out of place CROSS CORRELATION """

function crossCorrelation( ::ZNCC, g::A{<:Real,N}, f::A{<:Real,N}; corrType=Float32) where {N}

    sizef  = size( f ); 
    sizeg  = size( g );    
    
    csize   = sizef .+ sizeg; 
	shifts  = div.( csize, 2 ) .+ div.( size(f) .- size(g), 2 ) .- 1; 

    corr    = zeros( corrType, csize ); 
    shifted = zeros( corrType, csize ); 
    padf    = zeros( Complex{corrType}, csize ); 
    padg    = zeros( Complex{corrType}, csize );
    plan    =  plan_fft!( padf ); 
    iplan   = plan_ifft!( padf ); 

    
    meang, stdg = meanStd( g )
    meanf = Statistics.mean( f ); 
    
    putWithinPadded!( padf, f, meanf, (0,0,0), (1,size(f,1),1,size(f,2),1,size(f,3)) ); 
    putWithinPadded!( padg, g, meang, (0,0,0), (1,size(g,1),1,size(g,2),1,size(g,3)) );

    sumf2 = zeros( corrType, sizef .+ 1 ); 
    integralArraySQ!( padf, sumf2 );  
   
    
    crossCorrelation!( ZNCC(), corr, shifted, shifts, padf, padg, plan, iplan, 
                               sumf2, stdg, sizef, sizeg );
    
    return corr; 
end
