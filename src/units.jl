import Base.*, Base.show

abstract type DistanceUnits end

struct DistanceUnit <: DistanceUnits
    value::Float64
    factor::Int8
    name::String
end

DistanceUnit( fact::T1, name::String ) where {T1<:Integer} = DistanceUnit( 1.0, Int8(fact), name )
DistanceUnit( val::T1, fact::T2, name::String ) where {T1<:Real, T2<:Integer} = DistanceUnit( Float64(val), Int8(fact), name )


# defining and exporting each distance unit name
prefixes = [ ("yotta",24), ("zetta",21), ("exa",18), ("peta",15), ("tera",12), ("giga",9), ("mega",6),
             ("kilo",3), ("hecto",2), ("deca",1), ("",0), ("deci",-1), ("centi",-2), ("milli",-3),
             ("micro",-6), ("nano",-9), ("pico",-12), ("femto",-15), ("atto",-18), ("zepto",-21), ("yocto",-24) ]

for ( prefix, exponent ) in prefixes
    definition = "const $(prefix)m = DistanceUnit( $(exponent), \"$(prefix)meter\" )";
    exportstatement = "export $(prefix)m";

    eval( Meta.parse( definition ) );
    eval( Meta.parse( exportstatement ) );
end

Base.:*( num::T, un::DistanceUnits ) where {T<:Real} = DistanceUnit( num, un.factor, un.name )

changeUnits( from::DistanceUnits, to::DistanceUnits ) = DistanceUnit( from.value * 10.0^(from.factor - to.factor), to.factor, to.name )

Base.show(out::IO, un::DistanceUnits ) where U = print(out, un.value, un.name*"s" )


abstract type TimeUnits end

struct TimeUnit <: TimeUnits
    value::Float64
    seconds::Float64
    name::String
end

TimeUnit( secs::T1, name::String ) where {T1<:Real} = TimeUnit( 1.0, Float64(secs), name )
TimeUnit( val::T1, fact::T2, name::String ) where {T1<:Real, T2<:Real} = TimeUnit( Float64(val), Float64(fact), name )


# defining and exporting each time unit name
names = [ ("day",86400), ("hour",3600), ("minute",60), ("second",1), ("millisecond",0.001) ]

for ( name, numseconds ) in names
    definition = "const $(name) = TimeUnit( $(numseconds), \"$(name)\" )";
    exportstatement = "export $(name)";

    eval( Meta.parse( definition ) );
    eval( Meta.parse( exportstatement ) );
end

Base.:*( num::T, un::TimeUnits ) where {T<:Real} = TimeUnit( num, un.seconds, un.name )

changeUnits( from::TimeUnits, to::TimeUnits ) = TimeUnit( from.value * from.seconds/to.seconds, to.seconds, to.name )

Base.show(out::IO, un::TimeUnits ) where U = print(out, un.value, un.name*"s" )
