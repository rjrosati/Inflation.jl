using Inflation
using SymPy

using Random
using Roots
using ForwardDiff

d = 10

# sample masses from the Marchenko-Pastur distribution
m = 1e-5
β = 1/2
a = m^2 * (1-sqrt(β))^2
b = m^2 * (1+sqrt(β))^2
function MarchenkoPasturPDF(x)
    return sqrt( (b-x)*(x-a)) / (2*pi*x*β*m^2)
end
support = (a,b)

function rejection_sample(pdf,num,support=(-Inf,Inf))
    # for now assume pdf has exactly one critical point
    # will probably fail without finite support
    a = support[1]
    b = support[2]
    dpdf(x) = ForwardDiff.derivative(pdf,x)
    xmax = find_zero(dpdf,(a,b))
    pdf_max = pdf(xmax)
    samples = Array{Float64}(undef,num)
    i = 1
    while i <= num
        x = a+(b-a)*rand()
        s = rand()*pdf_max
        if s < pdf(x)
            samples[i] = x
            i += 1
        end
    end
    return samples
end

masses = rejection_sample(MarchenkoPasturPDF,d,support)

params = [symbols("m$i") for i in 1:d]
pvalues = masses

function G(Phi,p)
    d = length(Phi)
    g = Matrix{Float64}(I,d,d)
    return g
end


function V(Phi,p)
    m = p[1:end]
    v = sum( m .* Phi.^2 ) / 2
    return v
end

# pick an initial radius to give N0 e-folds
N0 = 100
# get a vector on the sphere, use gaussian trick
Phi0 = randn(d)
Phi0 ./= norm(Phi0)
Phi0 .*= sqrt(4*N0)

Pi0 = zeros(d)

funcs = inflation_setup(d,V,G,params=params)

sol = background_evolve(Phi0,Pi0,pvalues,funcs,verbose=true)

tsol = transport_perturbations(sol,pvalues,funcs,verbose=true)
println(tsol)

