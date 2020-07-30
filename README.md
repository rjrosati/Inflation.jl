# Inflation.jl

[![Build Status](https://travis-ci.com/rjrosati/Inflation.jl.svg?token=zMDX3GmCZbdBcf9JWMdp&branch=master)](https://travis-ci.com/rjrosati/Inflation.jl)
[![codecov](https://codecov.io/gh/rjrosati/Inflation.jl/branch/master/graph/badge.svg?token=FHiP2Mqh90)](https://codecov.io/gh/rjrosati/Inflation.jl)

# Inflation.jl
A Julia package for numerical evaluation of cosmic inflation models. Perturbations are evolved with the transport method.
![logo](https://github.com/rjrosati/Inflation.jl/raw/master/inflationjl.png "Logo")

## an example
```julia
using Inflation
using SymPy

# this is a quadratic inflation model with 10 fields and randomly selected masses

# set the field space dimension
d = 10

# set the model parameters, and their values
params = [symbols("m$i") for i in 1:d]
pvalues = rand(d)*1e-9

function G(Phi,p)
    d = length(Phi)
    g = Matrix{eltype(Phi)}(I,d,d)
    return g
end

function V(Phi,p)
    m = p
    return sum([ m[i]*Phi[i]^2)/2 for i in 1:d])
end

# set the initial conditions
Pi0 = zeros(d) # zero initial velocity

# pick a position randomly over the sphere, with radius to give ~N0 e-folds
N0 = 100
Phi0 = randn(d)
Phi0 ./= norm(Phi0)
Phi0 .*= sqrt(4*N0)

# symbolically construct the equations of motion
# because this potential is sum-separable, the horizon-crossing approximation will automatically be calculated
funcs = inflation_setup(d,V,G,params)

# solve the background equations of motion
sol = background_evolve(Phi0,Pi0,pvalues,funcs,verbose=true)

# solve the 2-pt correlation function equations of motion
# by default, use 7 k-values centered at a pivot scale of 0.002 Mpc^-1
tsol = transport_perturbations(sol,pvalues,funcs,verbose=true)

println(tsol)
```

see the `examples` directory for more complicated potentials and metrics, how to scan parameter space, evaluate several simulations in parallel, output and analyze data.

At the moment, `Inflation.jl` can only solve the 2-pt correlation function equations of motion. Solving for higher-point correlation functions is possible with [PyTransport](https://github.com/jronayne/PyTransport)/[CppTransport](https://github.com/ds283/CppTransport).
