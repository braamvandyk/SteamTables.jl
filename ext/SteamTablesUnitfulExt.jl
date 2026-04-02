module SteamTablesUnitfulExt

using Unitful
using Unitful: @u_str
using SteamTables
using SteamTables.PrecompileTools
import SteamTables: mysignif

struct UnitsError <: Exception
    var
    message::String
end


Base.showerror(io::IO, e::UnitsError) = print(io, "UnitsError with ",e.var, ": ", e.message)


function pt_unwrap_units(P::Quantity,T::Quantity)
    try
        _P = 1.0*uconvert(u"MPa", P)
        _T = 1.0*uconvert(u"K", T)
        return _P.val,_T.val
    catch
        throw(UnitsError((P,T), "Invalid input units."))
    end
end


function ph_unwrap_units(P::Quantity,h::Quantity)
    try
        _P = 1.0*uconvert(u"MPa", P)
        _h = 1.0*uconvert(u"kJ/kg", h)
        return _P.val,_h.val
    catch
        throw(UnitsError((P,h), "Invalid input units."))
    end
end


function ps_unwrap_units(P::Quantity,s::Quantity)
    try
        _P = 1.0*uconvert(u"MPa", P)
        _s = 1.0*uconvert(u"kJ/kg/K", s)
        return _P.val,_s.val
    catch
        throw(UnitsError((P,s), "Invalid input units."))
    end
end


macro func_def_with_unit(f,unit)
    f_Ps = Symbol(f, :_Ps)
    f_Ph = Symbol(f, :_Ph)
    ff = Symbol(f)
    quote
        function SteamTables.$f(P::Q1, T::Q2) where Q1 <: Quantity where Q2 <: Quantity
            _p, _T = pt_unwrap_units(P, T)
            return SteamTables.$f(_p, _T) * Unitful.@u_str($unit)
        end

        function SteamTables.$f_Ph(P::Q1, h::Q2) where Q1 <: Quantity where Q2 <: Quantity
            _P, _h = ph_unwrap_units(P, h)
            return SteamTables.$f_Ph(_P, _h) * Unitful.@u_str($unit)
        end


        function SteamTables.$f_Ps(P::Q1, s::Q2) where Q1 <: Quantity where Q2 <: Quantity
            _P, _s = ps_unwrap_units(P, s)
            return SteamTables.$f_Ps(_P, _s) * Unitful.@u_str($unit)
        end
    end |> esc
end


macro func_defT_with_unit(f,unit)
    quote 
        function SteamTables.$f(T::Q) where Q <: Quantity
            try
                _T = 1.0*uconvert(u"K", T).val
                return SteamTables.$f(_T)*Unitful.@u_str($unit)
            catch
                throw(UnitsError(T, "Invalid input units."))
            end
        end
    end |> esc
end


@func_def_with_unit SpecificF "kJ/kg"
@func_def_with_unit SpecificG "kJ/kg"
@func_def_with_unit SpecificU "kJ/kg"
@func_def_with_unit SpecificH "kJ/kg"
@func_def_with_unit SpecificS "kJ/kg/K"
@func_def_with_unit SpecificV "m^3/kg"
@func_def_with_unit SpecificCP "kJ/kg/K"
@func_def_with_unit SpecificCV "kJ/kg/K"
@func_def_with_unit SpeedOfSound "m/s"


function SteamTables.Tsat(P::Q) where Q <: Quantity
    try
        _P = 1.0*uconvert(u"MPa", P)
        return Tsat(_P.val)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
end


function SteamTables.Temperature_Ph(P::Q1, h::Q2) where Q1 <: Quantity where Q2 <: Quantity
    _P,_h = ph_unwrap_units(P, h)
    return Temperature_Ph(_P, _h)*u"K"
end


function SteamTables.Temperature_Ps(P::Q1, s::Q2) where Q1 <: Quantity where Q2 <: Quantity
    _P,_s = ps_unwrap_units(P, s)
    return Temperature_Ps(_P, _s)*u"K"
end


@func_defT_with_unit Psat "MPa"
@func_defT_with_unit SatDensL "kg/m^3"
@func_defT_with_unit SatDensV "kg/m^3"
@func_defT_with_unit SatHL "kJ/kg"
@func_defT_with_unit SatHV "kJ/kg"
@func_defT_with_unit SatSL "kJ/kg/K"
@func_defT_with_unit SatSV "kJ/kg/K"
@func_defT_with_unit DeltaHvap "kJ/kg"


@inline function SteamTables.mysignif(x::Q, n::Int) where Q <: Quantity
    z = x.val
    y = Quantity(round(z, sigdigits = n), unit(x))
    return y
end


if isdefined(Base,:get_extension)
    @setup_workload begin
        include("compilefile_Unitful.jl")

        @compile_workload begin
            let
                dummy = runprecompworkload_unitful()
                dummy = nothing
            end
        end
    end
end

end #module