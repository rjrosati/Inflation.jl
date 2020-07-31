using Inflation
using SymPy

# this is a helix-like potential, published in https://arxiv.org/abs/1905.07495

d = 3
function G(Phi,p)
    return [1 0 0 ; 0 1 0; 0 0 1]
end

params = [symbols(x) for x in split("z0 r del a f s l zf tf")]
pvalues = [0.397096578518034, 0.7, 2.0, 0.003, 0.0004, 0.001,1.5e-3,4.8,0.1]

function V(Phi,p)
    z0,R,Δ,A,f,σ,l,zend,tanf = p
    x,y,z = Phi
    v = l^4*(exp(z/R) + Δ*(1 - exp(-((A*cos(z/f)-x)^2 + (A*sin(z/f)-y)^2)/(2*σ^2))))
    return v
end

z0,R,Δ,A,f,σ,l,zend,tanf = pvalues
Phi0 = [A*cos(z0/f), A*sin(z0/f), z0]
c = Base.atan(-f^2+6*(A^2+f^2)*R^2,2*f*R)
b = A*f*sqrt(4*f^2*R^2+(f^2-6*(A^2+f^2)*R^2)^2)*σ^2 / ( (A^2+f^2)*R*(-f^2+6*(A^2+f^2)*R^2)*Δ)

zprime = -1/R/(1+A^2/f^2) #asymptotic z velocity
theta = z0/f+c
thetaprime = zprime
dr = b*exp(z0/R)
drprime = b/R*exp(z0/R)*zprime
xprime = -A/f*sin(z0/f)*zprime + drprime*cos(theta) - thetaprime*dr*sin(theta)
yprime =  A/f*cos(z0/f)*zprime + drprime*sin(theta) + thetaprime*dr*cos(theta)
Pi0 = [xprime,yprime,zprime]
println("Expected steady-state parameters")
eh0 = f^2 / (2*(A^2+f^2)*R^2)
println("ϵH = $eh0")
println("Ω  = $(9*4*A^2*f^2*R^2 / ((A^2+f^2)*(f^2-6*(A^2+f^2)*R^2)^2))")
println("ϵV = $(f^2*((A^2+f^2)/R^2 + 4*A^2*f^2/(f^2 -6*(A^2+f^2)*R^2)^2)/ (2*(A^2+f^2)^2))")
println("ηH = 0")
println("ns = $(1 - 2*eh0)")
println("r  = $(16*eh0)*cs")

funcs = inflation_setup(d,V,G,params)

println("Starting background evolution")

options = BackgroundOptions(
                            max_efolds=75.0
                           )
sol = background_evolve(Phi0,Pi0,pvalues,funcs,verbose=true,options=options)


tsol = transport_perturbations(sol,pvalues,funcs,verbose=true)
output_data("helix.jld2",funcs,pvalues,sol,tsol)


