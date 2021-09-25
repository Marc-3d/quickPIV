#=
	we make heavy use of custom type aliases. Have a look at the first lines of quickPIV.jl, 
	to see the type abreviations used throughout the file "III", "I", "F"
=#

""" PIV PARAMETER OBJECT """

mutable struct PIVParameters

         corr::CORRTYPES # cross-correlation algorithm
    interSize::III       # interrogation area/volume size in pixels/voxels
 searchMargin::III       # search margin in pixels/voxels
      overlap::III       # overlap in pixels/voxels between adjacent IA/IV
        mpass::I         # multi-pass depth
        width::I         # exclusion radius around 1st peak for finding 2nd peak
         peak::String    # one of "", "gaussian", "centroid"
     sigNoise::String    # one of "", "peak2peak"
      filtFun::Function  # Filtering function to detect IA/IV corresponding to background
    threshold::F         # skip cross correlation if filtFun( interArea ) < threshold
        units::Tuple{DU,DU,DU,TU} # spatial and temporal units of the analyzed data
         args::Tuple{III,III,III,I,I,String,String,Function,F}
end

function parseCorr( corr::String )
    if     lowercase(corr) == "zncc"
		return   ZNCC();
	elseif lowercase(corr) == "fft"
		return    FFT();
	elseif lowercase(corr) == "nsqecc" || lowercase(corr) == "diff"
		return NSQECC();
	else
		return NSQECC();
	end
end
toIII( input::I   ) = ( input, input, input )
toIII( input::II  ) = ( input[1], input[2], 0 )
toIII( input::III ) = input

# convenient construction that takes care of filling up "args"
function PIVParameters( corr::CORRTYPES, interSize::UnI, searchMargin::UnI, ovp::UnI, mpass::I,
                        width::I, peak::String, SN::String, ffn::Function, th, uns::UNTS )

	th_f = th + 0.0; # converting to float
	args = (toIII(interSize), toIII(searchMargin), toIII(ovp), mpass, width, peak, SN, ffn, th_f)

	return PIVParameters( corr, args..., uns, args )
end

# constructor with default values
function setPIVParameters(;corr = "nsqecc", interSize = 32, searchMargin = 0, overlap = 0,
                           mpass = 1, width = 2, peak = "gaussian", sigNoise = "", 
                           filtFun = (x)->maxval(x), threshold = -1.0, units = (1m,1m,1m,1second))

    return PIVParameters( parseCorr(corr), interSize, searchMargin, overlap, mpass,
                          width, peak, sigNoise, filtFun, threshold, units )
end

function checkPIVParameters( dsize::NTuple{N,I}, params ) where {N}

    !( N !== 2  ||  N !== 3 ) && error("Your input data isn't 2D nor 3D.")

	any( params.interSize[1:N] .> dsize ) && error("interSize is bigger than the data")
    any( params.interSize .< 2 ) && error("interSize = 1 makes little sense.")
    any( params.overlap .>= params.interSize ) && error("Overlap must be smaller than interSize")
end



""" SYNTHETIC PARTICLES PARAMETER OBJECT """

mutable struct synthParameters
      n::I    # total number of particles, ignored when dens > 0
   dens::I    # particle density ( nºparticles/IA )
      d::I    # depth. If depth > 1, generates a volume
      w::I    # width of generated image/volume
      h::I    # height of generated image/volume
     i0::F    # max intensity of particles
     dt::F    # standard deviation of gaussian around each particles
     th::F    # standard deviation of gaussian around lightsheet
    err::F    # determines rendering radius around each pixel
      z::I    # light-sheet position. Default z = 0
  noise::I    # varible for controlling noise. Not really implemented.
    rad::I    # overwrites the rendering radius computed from err.
   mode::Bool # switches between ( x, y, z ) or ( y, x, z ) coordinates
   args::Tuple{I,I,I,I,I,F,F,F,F,I,I,I,Bool}
end

# convenience constructor
function synthParameters( n::I, dens::I, d::I, w::I, h::I, i0::F, dt::F, th::F, err::F, z::I, 
                          noise::I, rad::I, mode::Bool )

	args = (n, dens, d, w, h, i0, dt, th, err, z, noise, rad, mode)
	return synthParameters( args..., args )
end

# constructor with default values
function setSynthParameters(; n=0::I, dens=0::I, d=1::I, w=300::I, h=300::I,
                              z=0::I, i0=255.0::F, dt=3.0::F, err=0.1::F,
                              th=10.0::F, noise=0::I, rad=0::I, mode="whd"::S )

    return synthParameters( n, dens, d, w, h, i0, dt, th, err, z, noise, rad, mode=="whd" )
end



""" TRANSFORMATION PARAMETERS """

mutable struct transformParameters{N,M}
      kind::Symbol      # currently ony :translation is supported
     means::NTuple{N,F} # mean value of transform in x, y [ and z ] if N = 3
      vars::NTuple{N,F} # variance component in covariance matrix
 cov_ratio::NTuple{M,F} # covariance components in covariance matrix. -1 < 1
      args::Tuple{Symbol,NTuple{N,F},NTuple{N,F},NTuple{M,F}}
end

TP2D = transformParameters{2,1}
TP3D = transformParameters{3,3}

# convenience constructor
function transformParameters( kind::Symbol, means::NTuple{N,F}, vars::NTuple{N,F},
                              covs::NTuple{M,F} ) where {N,M}

	args = ( kind, means, vars, covs )
	return transformParameters{N,M}( args..., args )
end

# constructor with default values
function setTransformParameters(; kind=:translation, means=(1.0, 1.0), vars=(eps(Float64),
								  eps(Float64)), cov_ratio=(0.0,) )

    return transformParameters( kind, means, vars, cov_ratio );
end

to2DTransform( t2D::TP2D ) = t2D

function to2DTransform( t3D::TP3D )
    means     = (  t3D.means[1], t3D.means[2] );
    vars      = (  t3D.vars[1] ,  t3D.vars[2] );
    cov_ratio = (  t3D.cov_ratio[1], );

    return transformParameters{2,1}( t3D.kind, means, vars, cov_ratio );
end

function to3DTransform( t2D::TP2D )
    means     = (   t2D.means[1]  ,   t2D.means[2]  , (t2D.means[1] + t2D.means[2])/2 );
    vars      = (    t2D.vars[1]  ,    t2D.vars[2]  , ( t2D.vars[1] +  t2D.vars[2])/2 );
    cov_ratio = ( t2D.cov_ratio[1], t2D.cov_ratio[1], t2D.cov_ratio[1] );

    return transformParameters{3,3}( t2D.kind, means, vars, cov_ratio );
end

to3DTransform( t3D::TP3D ) = t3D



""" EVALUATION PARAMETER OBJECT """

mutable struct metaParameters
 metaloop::I      # determines the number of iterations between "min" and "max" values
 variable::Symbol # specifies what variable will be changed each iteration.
      min::F      # minimum value of changing variable
      max::F      # maximum value of changing variable
  repeats::I      # number of measurements for computing bias and error
     args::Tuple{I,Symbol,F,F,I}
end

function metaParameters( metaloop::I, variable::Symbol, min::F, max::F, repeats::I )

	args = (metaloop, variable, min, max, repeats )
	return metaParameters( args..., args )
end

# short type aliases, to make function signatures shorter
PP   = PIVParameters;
SP   = synthParameters;
TP   = Union{TP2D,TP3D}
MP   = metaParameters;


"""
    Creates a set of parameters for directing PIV evaluation. Accepted parameters are:
        variable (  Symbol ): a symbol indicating what variable should be evaluated. Can be any parameters
                              in synthParameters, PIVParameters or transformParameters.
        metaloop (  Int64  ): number of iterations
            min  ( Float64 ): min value of changing variable
            max  ( Float64 ): max value of changing variable
        repeats  (  Int64  ): number of repetitions
"""
function setMetaParameters( ; metaloop=0::I, var=Symbol("")::Symbol, min=0.0::F,  max=0.0::F, repeats=1::I )

    return metaParameters( metaloop, var, min,  max, repeats )
end


# function to update parameters, and also update the .args tuple, which was created to pass
# parameters by tuple deconstruction "parameter.args..."

function updateParameters!( params, field, newVal )
    hasfield = false;
    fieldidx =     0;
    for e in fieldnames( typeof(params) )
        fieldidx += 1;
        if ( e == field )
            hasfield = true;
            setproperty!( params, e, newVal )
            break
        end
    end

    ( !hasfield ) && return nothing;

    if typeof(params) == PIVParameters
        fieldidx -= 1 # ignore corr parameter
    end
    newargs = [ arg for arg in params.args ]
    newargs[ fieldidx ] = newVal;
    setproperty!( params, :args, Tuple(newargs) )
end
