using Inflation
using SymPy

# this is a 2-field non-minimally coupled model, published in https://arxiv.org/pdf/1304.0363.pdf

d = 2

function G(Phi,p)
    f = F(Phi,p)
    id = Matrix{Int64}(I,d,d)
    gf = [diff(f,phi) for phi in Phi]
    g = 1/(2*f) .* (id .+ 3*gf'.*gf./f)
    return g
end

xis = [symbols("xi$i",positive=true) for i in 1:d]
lambdas = [symbols("l$i",positive=true) for i in 1:d]
params = [xis; lambdas; symbols("g",positive=true) ]
pvalues = [ 1e3, 1.2*1e3, 1e-2, 0.75*1e-2, 1e-2]

function F(Phi,p)
    xi = p[1:d]
    return 1/2*(1 + xi' * (Phi.^2))
end

function V(Phi,p)
    λs = p[d+1:end-1]
    g = p[end]
    jordanV = g/2*prod(Phi.^2) + λs'*(Phi.^4)/4
    f = F(Phi,p)
    v = jordanV / (4*f^2)
    return v
end

funcs = inflation_setup(d,V,G,params)


ξϕ,ξχ = pvalues[1:2]
λϕ,λχ = pvalues[3:4]
r0 = 10*max(1/sqrt(ξϕ),1/sqrt(ξχ))

# scan the number of e-folds, as a function of angle in the ϕ, χ plane
ts = range(0,2*pi,length=100)
Nends = []
options = BackgroundOptions(
                            # these trajectories can have transients with very high ϵₕ
                            #allow 1 efold for them to clear
                            min_efolds = 1.0,
                            max_efolds = 805.0
                           )
for t0 in ts
    println("r0: $r0")
    println("t0: $t0")
    Phi0 = r0*[cos(t0),sin(t0)]
    Pi0 = zeros(d) # works well enough


    sol = background_evolve(Phi0,Pi0,pvalues,funcs,verbose=true,options=options)
    push!(Nends,sol.t[end])
end


using Plots
using LaTeXStrings
plot(ts,Nends,xlabel=L"\theta",ylabel=L"N_e")

