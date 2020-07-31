# this potential is an α-attractor model with a quadratic potential
# this run script is set to run on a cluster with 24 total physical CPU cores
# the simulations are scanned over the model, for various values of α, N0, and Nf
# each simulation, if it successfully computes perturbations is saved in a data folder

const NUMCPU = 24
Nfs = [2,5,10]
runs = 100
alphas = [1//6, 1//60]
N0s = 75:25:400

datafolder = "/data/ANflation/"
subfolder = "batch1"

if !isdir(datafolder)
    println("Data folder doesn't exist!")
elseif !isdir(joinpath(datafolder,subfolder))
    println("Creating $(joinpath(datafolder,subfolder)) ...")
    mkdir(joinpath(datafolder,subfolder))
end


using Distributed
if length(workers()) == 1
    addprocs(NUMCPU - 1)
end

@everywhere using Inflation
@everywhere using SymPy



@everywhere begin
    using Random
    using Roots
    using ForwardDiff
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
end

for d in Nfs
    function G(Phi,p)
        d = length(Phi)
        α = p[1]
        Phi2 = sum([f^2 for f in Phi])
        g = Diagonal(ones(Sym,d) .* 6*α / (1 - Phi2)^2)
        return g
    end

    function V(Phi,p)
        m = p[2:end]
        α = p[1]
        v = 0
        for i in 1:d
            v += α/2.0*(m[i]*Phi[i]^2)
        end
        return v
    end

    runs_to_run = []
    for N0 in N0s, alpha in alphas, run in 1:runs
        fname = "N$(N0)_a$(Int(1/alpha))_Nf$(d)_r$run.jld2"
        fullfname = joinpath(datafolder,subfolder,fname)
        if !isfile(fullfname)
            push!(runs_to_run,(N0,alpha,fullfname))
        end
    end

    if length(runs_to_run) == 0
        println("All runs done for d=$d")
        continue
    end
    println("Running $(length(runs_to_run)) runs at d=$d")
    # all params should be lowercase
    mparams = [symbols("m$i",real=true,positive=true) for i in 1:d]
    params = [symbols("α",real=true,positive=true),mparams...]

    quotes = inflation_setup(d,V,G,params=params,attempt_HCA=false,return_quotes=true)
    @everywhere funcs = Inflation.make_funcs($quotes)
    @everywhere runs_to_run = $runs_to_run
    @sync @distributed for (N0,alpha,fullfname) in runs_to_run
        println(fullfname)
        while !isfile(fullfname)
            mass_eigs = rejection_sample(MarchenkoPasturPDF,d,support)
            pvalues = [alpha,mass_eigs...]

            r0 =sqrt( (2*N0 + 3*sqrt(alpha) - 3*alpha) / (2*N0 + 3*sqrt(alpha)))
            ## get a vector on the sphere, use gaussian trick
            # N0 prior should be approximately true
            Phi0 = randn(d)
            Phi0 ./= norm(Phi0)
            Phi0 .*= r0

            Pi0 = zeros(d)
            rNe = 3*alpha/2*(1/(1-Phi0'*Phi0)-1/sqrt(alpha))
            println("estimated radial efolds: $rNe")
            println("total energy: $(V(Phi0,pvalues))")

            #analytics
            if d==2
                θ0 = Base.atan(Phi0[2],Phi0[1])
                R = pvalues[end]/pvalues[end-1]
                a = pvalues[1]
                println("R: $R, α: $a, θ0: $θ0")
                Nang = -1/(54*a)*( (18*a+4)*log(0.5*(-(R-1)*cos(2*θ0)+R+1)) - 2*(R+2) + 8*(-(R^2-1)*cos(2*θ0)+R^2+R+1) / (-(R-1)*cos(2*θ0)+R+1)^2 )
                println("Nang: $Nang")
                Neff = Nexit_to_end - Nang
                if Neff < 0
                    println("Exit horizon during angular motion, no analytic result.")
                else
                    anal_ns = 1 - 2/Neff
                    anal_r  = 12*a/Neff^2
                    println("Expected ns: $anal_ns")
                    println("Expected r: $anal_r")
                end
            else
                println("d>2, no analytic result")
            end

            sol = background_evolve(Phi0,Pi0,pvalues,funcs,verbose=true)
            Nend = sol.t[end]
            Nearly = Nend - Nexit_to_end - Nsubhorizon - Nkrange
            if Nearly > 0
                tsol = transport_perturbations(sol,pvalues,funcs,verbose=true)
                output_data(fullfname,funcs,pvalues,sol,tsol)
            else
                println("Not enough inflation! retrying...")
                continue
            end
        end
    end
end
rmprocs(NUMCPU-1)
