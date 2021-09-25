VTK2Julia = Dict( "unsigned_char"  =>  UInt8, "char"  =>  Int8,
                  "unsigned_short" => UInt16, "short" => Int16,
                  "unsigned_int"   => UInt32, "int"   => Int32,
                  "unsigned_long"  => UInt64, "long"  => Int64,
                  "float" => Float32, "double" => Float64     )

Julia2VTK = Dict( UInt8   => "unsigned_char" , Int8  =>  "char",
                  UInt16  => "unsigned_short", Int16 => "short",
                  UInt32  => "unsigned_int"  , Int32 =>   "int",
                  UInt64  => "unsigned_long" , Int64 =>  "long",
                  Float32 => "float", Float64 => "double"      )

function parseFilename( filename, path, format )
    !( occursin( '/', filename ) ) && ( filename = path*filename );
    ( length(filename) < 4 ) && ( filename = filename*format )
    !( occursin( '.', filename[end-4:end] ) ) && ( filename = filename*format )
    return filename;
end

""" Loads 3D tiff stacks, 3D vtk volumes and any file openable with FileIO.load() """

function PIVload( filename::String; typ=nothing, sample=1, path=pwd() )

    fn = parseFilename( filename, path, "" )
    isfile(fn) || error("$fn can not be found")

    if any( occursin.( (".tif",".TIF",".TIFF",".tiff"), fn ) )
        data = LIBTIFF.tiffread( fn, typ=typ, sample=sample )

    elseif occursin( ".vtk", fn )
        data = loadVTKVolume( fn, typ=typ )
    else
        data = FileIO.load(fn)
    end

    return data
end

""" Storing 3D Vector fields to visualize in Paraview """

function vectorFieldToVTK( filename::String, components::Array{T,3}...;
                           path="", mode="w" ) where {T<:Real}

    data_type = get( Julia2VTK, T, "" );
    ( data_type == "" ) && ( error("Unrecognized data type") );

    dims    = length(components);
    U, V, W = components[ 1 ], components[ 2 ], components[ 3 ];
    h, w, d =   size( U );
    npoints = length( U );

    fn = parseFilename( filename, path, ".vtk" );

    io = open( fn, mode )

        println( io, "# vtk DataFile Version 2.0"  )
        println( io, "PIV3D.jl vector field"  )
        println( io, "ASCII"  )
        println( io, "DATASET STRUCTURED_GRID"  )
        println( io, "DIMENSIONS $w $h $d"  )
        println( io, "POINTS $npoints int"  )

        # from line 7 to line 7 + length(U) we store grid coordinates
		for z in 1:d
			for y in 1:h
				for x in 1:w
            		println( io, string( x, " ", y, " ", z ) )
				end
			end
        end

        println( io, "POINT_DATA $npoints"  )
        println( io, "VECTORS directions $data_type"  )

        # from line 7 + length(U) + 2 to line 7 + 2 + 2*length(U) we store vector coordinates
		for z in 1:d
			for y in 1:h
				for x in 1:w
            		println( io, string( V[y,x,z], " ", U[y,x,z], " ", W[y,x,z] ) )
				end
			end
        end

    close(io)
    println( "3D vtk vector field saved.")
end

""" Storing 3D volumes to visualize in Paraview """

function volumeToVTK( filename::String, volume::Array{T,3}; 
					  path="", mode="w", off=(0,0,0) ) where {T<:Real}

    data_type = get( Julia2VTK, T, "" );
    ( data_type == "" ) && ( error("Unrecognized data type") );

    h, w, d = size( volume )

    fn = parseFilename( filename, path, ".vtk" );
    io = open( fn, mode )

        println( io, "# vtk DataFile Version 2.0"  );
        println( io, "PIV3D.jl volume"  );
        println( io, "ASCII"  );
        println( io, "DATASET STRUCTURED_POINTS"  );
        println( io, "DIMENSIONS $w $h $d"  );
        println( io, "ORIGIN $(off[1]) $(off[2]) $(off[3])"  );
        println( io, "SPACING 1 1 1"  );
        println( io, "POINT_DATA $(length(volume))"  );
        println( io, "SCALARS intensities $data_type"  );
        println( io, "LOOKUP_TABLE default" );

        # from line 11 to line 11 + length(volume) we store each voxel intensity
		for z in 1:d
			for y in 1:h
				for x in 1:w
		    		println( io, string( volume[y,x,z] ) )
				end
			end
		end

    close(io)
    println( "VTK volume saved.")
end

""" Storing pseudo-trajectories to visualize in Paraview """

function trajectoriesToVTK( filename::String,
                            V::Array{T,2}, L::Array{T,2}, D::Array{T,2};
                            path="", mode="w" ) where {T<:Real}

    data_type = get( Julia2VTK, T, "" );
    ( data_type ==  "" ) && ( error("Unrecognized data type") );

    times, num_points = size( V )

    fn = parseFilename( filename, path, ".vtk" );

    io = open( fn, mode )

        println( io, "# vtk DataFile Version 2.0" )
        println( io, "PIV3D Trajectories" )
        println( io, "ASCII" )
        println( io, "DATASET POLYDATA" )
        println( io, "POINTS $(num_points*times) $(data_type)")
        for p in 1:num_points, t in 1:times
            println( io, string( "$(L[t,p]) $(V[t,p]) $(D[t,p])" ) )
        end

        println( io, "LINES $(num_points) $((times+1)*num_points)" )
        println( io, "" )
        idx = 0;
        for p in 1:num_points
            print( io, "$(times) " )
            for t in 1:times
                print( io, "$idx " )
                idx += 1
            end
            println( io, "" )
        end

        println( io, "POINT_DATA $(num_points*times)" )
        println( io, "SCALARS index int 1" )
        println( io, "LOOKUP_TABLE default" )
        println( io, "" );
        for p in 1:num_points, t in 1:times
            print( io, "$t " )
        end

    close(io)

    println( "trajectories vtk saved.")
end


""" loading VTK volumes stored with the PIV3D.volumeToVTK(...) function. """

function loadVTKVolume( filename::String; typ=nothing, path=pwd() )

    fn = parseFilename( filename, path, "" )

    lines = readlines( fn );
    if ( lines[2] !== "PIV3D.jl volume" )
		error("The vtk volume was not stored with PIV3D, and cannot open it" )
	end

    w,  h,  d = parse.( Int32 , split( lines[5], " " )[2:4] );
    data_type = ( typ == nothing ) ? get( VTK2Julia, dtype, Bool ) : typ;
    volume    = Array{ data_type, 3 }( undef, h, w, d )

    start_data, offset = 11, 0;
    for z in 1:d, r in 1:h, c in 1:w
		volume[r,c,z] = parse( data_type, lines[ start_data + offset ] )
		offset += 1;
    end
    return volume;
end

function VTKVectorFieldSize( filename::String; path=pwd() )

    fn = parseFilename( filename, path, ".vtk" );
    isfile(fn) || error("$fn can not be found")

    lines = readlines(fn);
    return parse.( Int32, split(lines[5]," ")[2:4] )
end


""" loading VTK vector fields stored with the PIV.vectorFieldToVTK(...) function """

function loadVTKVectorField( filename::String; typ=nothing, path=pwd() )

    fn = parseFilename( filename, path, ".vtk" )

    lines = readlines(fn);
    if ( lines[2] !== "PIV3D.jl vector field" )
		error("The vtk vector field was not stored with PIV3D, and cannot open it" )
	end

    w, h, d = parse.( Int32, split(lines[5]," ")[2:4] )

    _, _, dtype = split( lines[7+w*h*d+1], " " );
    data_type = ( typ == nothing ) ? get( VTK2Julia, dtype, Bool ) : typ;

    U = Array{data_type, 3}(undef, h, w, d);
    V = Array{data_type, 3}(undef, h, w, d);
    W = Array{data_type, 3}(undef, h, w, d);

    vector_start, offset = 7 +  w*h*d + 2, 0;

	for z in 1:d, y in 1:h, x in 1:w
		_v, _u, _w = parse.( data_type, split( lines[vector_start+offset] ) )
		U[y, x, z] = _u
		V[y, x, z] = _v
		W[y, x, z] = _w
		offset += 1;
	end

    return U, V, W
end
