module Inflation

using SymPy
@info "Travis CI hack"
using DifferentialEquations
@info "Travis CI hack"
using ODEInterfaceDiffEq
@info "Travis CI hack"
#using DiffEqParamEstim
using LinearAlgebra
using CurveFit
using Serialization
using JLD2
using Parameters

export inflation_setup
export background_evolve
export transport_perturbations
export SetupOptions, BackgroundOptions, PerturbationOptions
export output_data


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
