module Inflation

using SymPy
using DifferentialEquations
using ODEInterfaceDiffEq
#using DiffEqParamEstim
using LinearAlgebra
using CurveFit
using Serialization
using JLD2
using Parameters
using RuntimeGeneratedFunctions

export inflation_setup
export background_evolve
export transport_perturbations
export SetupOptions, BackgroundOptions, PerturbationOptions
export output_data


SymPy.diff(x::T,::Sym) where {T <: Real} = 0
SymPy.diff(x::T,::Sym,::Sym) where {T <: Real} = 0

RuntimeGeneratedFunctions.init(@__MODULE__)

# needed as of SymPy 1.0.38
import SymPy.fn_map
fn_map["Mul"] = :__prod__

# change this whenever there's a breaking change in the cache functions
cache_version = 1.3
# unit conversion
Mpc_to_Lpl = 2.6245e-57

include("options.jl")
include("quoting.jl")
include("symbolic_setup.jl")
include("solvers.jl")
include("output.jl")

end # module
