module quickPIV

using FFTW, FileIO, LIBTIFF, Statistics

export PIV

# Shortened type names to keep function signatures short ( when possible one-liners )
include("units.jl")

F      = Union{Float64,Float32} # considering 64 and 32 bit systems
I      = Union{Int64,Int32}
A      = Array
S      = String
C      = Complex
II     = Tuple{I,I}
III    = Tuple{I,I,I}
UnI    = Union{I,II,III}
SySySy = Tuple{Symbol,Symbol,Symbol}
DU     = DistanceUnits
TU     = TimeUnits
UNTS   = Tuple{DU,DU,DU,TU}

Base.convert( ::Type{Union{Int32,Int64}}, a::AbstractFloat ) = round( typeof(0), a )

struct ZNCC   end
struct FFT    end
struct NSQECC end

CORRTYPES = Union{ZNCC,FFT,NSQECC}; 

include("parameters.jl")
include("pivIO.jl")
include("utils.jl")

include("PIV/subpixel.jl")
include("PIV/signalToNoise.jl")
include("PIV/ccrZNCC.jl")
include("PIV/ccrFFT.jl")
include("PIV/ccrNSQECC.jl")

include("syntheticParticles.jl")

include("postprocessing.jl")


function PIV( data1::Array{T,N}, data2::Array{T,N}, params::PIVParameters;
              checkParams=true, units=nothing ) where {T<:AbstractFloat,N}

	checkParams && checkPIVParameters( size(data1), params )

	return _PIV( data1, data2, params )
end

""" Calling 2D PIV analyses """

function _PIV( img1::A{T,2}, img2::A{T,2}, p::PIVParameters; units=nothing 
             ) where {T<:AbstractFloat}

	u, v, sn = PIV_2D( p.corr, img1, img2, p.args... );

	if units !== nothing
	    ratios = changeUnits.(  p.units[1:2], units[1:2] )./changeUnits( p.units, units[4]  )
	    for e in 1:length(u)
	        u[e] = u[e]*ratios[1]
	        v[e] = v[e]*ratios[2]
	    end
	end
	return u, v, sn;
end

""" Calling 3D PIV analyses """

function _PIV( vol1::A{T,3}, vol2::A{T,3}, p::PIVParameters; units=nothing 
             ) where {T<:AbstractFloat}

	u, v, w, sn = PIV_3D( p.corr, vol1, vol2, p.args... );

	if units !== nothing
	    ratios = changeUnits.( p.units, units[1:3] )./changeUnits( p.units, units[4] )
	    for e in 1:length(u)
	        u[e] = u[e]*ratios[1]
	        v[e] = v[e]*ratios[2]
	        w[e] = w[e]*ratios[3]
	    end
	end

	return u, v, w, sn;
end

end
