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

export inflation_setup
export background_evolve
export transport_evolve
export SetupOptions
export SolverOptions
export output_data


# change this whenever there's a breaking change in the cache functions
cache_version = 1.2
# unit conversion
Mpc_to_Lpl = 2.6245e-57

include("symbolic_setup.jl")
include("quoting.jl")
include("options.jl")
include("solvers.jl")
include("output.jl")

end # module
