function background_evolve(Phi0,Pi0,pvals,funcs;verbose=false,
    options::BackgroundOptions = BackgroundOptions())
    Eh = funcs["Eh"]
    accel = funcs["Pi_eom"]
    V = funcs["V"]
    Hub = funcs["H"]
    d = length(Phi0)
    function eom(du,u,p,t)
        # name fields
        Phi  = u[1:d]
        Pi   = u[d+1:2d]
        par  = [Phi;Pi;pvals]
        #e    = Eh(par...)
        # eom
        du[1:d] = Pi
        du[d+1:2*d] = accel(par...)
    end
    function end_inflation(vars,t,integrator)
        if t > options.min_efolds
            e = Eh(vars...,pvals...)
            return options.end_epsilon-e # trigger when this hits zero from above
        else
            return 1.0
        end
    end
    e0 = Eh([Phi0;Pi0;pvals]...)
    if e0 >= 1.0 && iszero(options.min_efolds)
        #println("Initial ϵH = $e0")
        error("Inflation cannot start, initial eH = $e0 > 1")
    end
    H0 = Hub([Phi0;Pi0;pvals]...)
    if verbose
        vars = [Phi0;Pi0;pvals]
        ev0 = funcs["Ev"](vars...)
        et0 = funcs["Et"](vars...)
        o0 = funcs["Om"](vars...)
        println("Initial H: $H0")
        println("Initial ϵV: $ev0")
        println("Initial ϵH: $e0")
        println("Initial η: $et0")
        println("Initial Ω: $o0")
    end

    func = ODEFunction(eom)#;jac=eom_jac)
    SolT = eltype(Phi0)
    prob = ODEProblem(func,SolT[Phi0;Pi0],(zero(SolT),SolT(options.max_efolds)))
    cb = ContinuousCallback(end_inflation,terminate!)
    if verbose
        println("Solving background eom")
    end
    sol = DifferentialEquations.solve(prob,options.background_solver,callback=cb,
    reltol=options.background_rtol,abstol=options.background_atol,
    maxiters=options.max_steps,
    save_everystep=options.save_everystep,saveat=options.saveat)# automatically choose algorithm?
    if verbose
        println("Done. Nend = $(sol.t[end])")
    end
    Nend = sol.t[end]
    if Nend > options.Nexit_to_end && haskey(funcs,"ns_hca") && funcs["ns_hca"] != Nothing
        vars = [sol(Nend-options.Nexit_to_end);pvals]
        println("Estimated Pζ (HCA):\n $(funcs["Pz_hca"](vars...))")
        println("Estimated ns (HCA):\n $(funcs["ns_hca"](vars...))")
    end

    return sol
end

function geodesic(Phi0,Pi0,pvals,funcs;verbose=false,max_efolds=max_efolds)
    accel = funcs["Geo_eom"]
    function eom(du,u,p,t)
        # name fields
        Phi  = u[1:d]
        Pi   = u[d+1:2d]
        par  = [Phi;Pi;pvals]
        # eom
        du[1:d] = Pi
        du[d+1:2*d] = accel(par...)
    end

    func = ODEFunction(eom)
    prob = ODEProblem(func,Float64[Phi0;Pi0],(0.0,max_efolds))
    if verbose
        println("Solving background eom")
    end
    sol = DifferentialEquations.solve(prob,background_solver,reltol=background_rtol,abstol=background_atol,maxiters=max_steps,save_everystep=false,saveat=0.001)# automatically choose algorithm?
    if verbose
        println("Done. Nend = $(sol.t[end])")
    end

    return sol

end

function get_ks(sol,funcs,pvals,options)
    d = length(sol.u[1])÷2
    Hub(N) = funcs["H"]([sol(N);pvals]...)
    Nend = sol.t[end]
    Npivot = Nend - options.Nexit_to_end
    loga_init = log(options.kpivot) + log(Mpc_to_Lpl) - log(Hub(Npivot)) - Npivot

    Nkrange = options.Nkrange
    Nstep = 2*Nkrange / options.desired_ksteps
    Ns = Npivot-Nkrange:Nstep:Npivot+Nkrange
    Hs = [Hub(n) for n in Ns]
    logks = Ns .+ loga_init .+ log.(Hs)
    return loga_init,logks,Ns
end

function curved_nullspace(g,v)
    d = length(v)
    proj(u,v) = u * (u'*g*v) / (u'*g*u)
    norm(v) = sqrt(v'*g*v)
    ginv = inv(g)
    v ./= norm(v)
    perpuu = ginv .- v'.*v
    ws = [ perpuu[:,i] ./ norm(perpuu[:,i]) for i in 2:d ]
    ortho_ws = Array{eltype(v),2}(undef,d,d)
    ortho_ws[:,1] = v
    for i in 2:d
        w = ws[i-1]
        for j in 1:i-1
            w .-= proj(ortho_ws[:,j],w)
        end
        w ./= norm(w)
        ortho_ws[:,i] = w
    end
    return ortho_ws[:,2:end]
end

function transport_perturbations(sol,pvals,funcs;verbose=false,options::PerturbationOptions =PerturbationOptions())
    d = length(sol.u[1])÷2
    Nend = sol.t[end]
    Nsubhorizon = options.Nsubhorizon
    Nexit_to_end = options.Nexit_to_end
    Nkrange = options.Nkrange
    Nearly = Nend - Nexit_to_end - Nsubhorizon - Nkrange
    if Nearly < 0
        error("Not enough inflation! Attempted to evolve perturbations, but would've started at Ne = $Nearly")
    end
    logainit,logks,Nexits = get_ks(sol,funcs,pvals,options)
    if verbose
        println("Solving perturbations, with $(length(logks)) k-values")
        println("Nexit_to_end: $Nexit_to_end")
    end
    H(N) = funcs["H"]([sol(N);pvals]...)
    Mab(N) = funcs["Mab"]([sol(N)[1:2*d];pvals]...)
    Eh(N) = funcs["Eh"]([sol(N)[1:2*d];pvals]...)
    Ev(N) = funcs["Ev"](sol(N)...,pvals...)
    Vab(N) = reshape(funcs["Vab"](sol(N)...,pvals...),d,d)
    Om(N) = funcs["Om"]([sol(N)[1:2*d];pvals]...)
    Om2(N) = funcs["Om2"]([sol(N)[1:2*d];pvals]...)
    ginv(N) = funcs["ginv"]([sol(N)[1:2*d];pvals]...)
    g(N) = funcs["g"]([sol(N)[1:2*d];pvals]...)
    #eom = funcs["Pert_eom"]
    uab(t,loga,logk) = reshape(funcs["Pert_uab"](t,loga,logk,[sol(t)[1:2*d];pvals]...),2*d,2*d)
    hp(t) = reshape(funcs["Pert_curv"]([sol(t)[1:2*d];pvals]...),d,d)
    #jac = funcs["Pert_jac"]
    #function Peom_jac(GJ,G,logk,t)
    #    GJ[:,:] = reshape(jac(t,logainit,logk,G...,[sol(t)[1:2*d];pvals]...),(2*d)^2,(2*d)^2)
    #    nothing
    #end
    function Gprime(dG,G,logk,t)
        @assert(1.0+1e-8 >= Eh(t) >= 0.0, "ϵH weird value: N=$t, ϵH=$(Eh(t))")
        Gπ = hp(t)
        G11 = @view(G[1:d,:])
        G12 = @view(G[d+1:end,:])
        #G21 = @view(G[:,1:d])
        #G22 = @view(G[:,d+1:end])
        dG .= uab(t,logainit,logk)*G .- vcat(Gπ*G11, Gπ*G12)# .+ hcat(G21*Gπ, G22*Gπ)
        #dG[:,:] = reshape(eom(t,logainit,logk,G...,[sol(t)[1:2*d];pvals]...),2*d,2*d)
        nothing
    end
    Nend = sol.t[end]
    Nhc = sol.t[end] - Nexit_to_end
    Gs = Function[]
    odefunc = ODEFunction(Gprime)#,jac=Peom_jac)
    SolT = eltype(sol(0.0))
    for (logk,Nexit) in zip(logks,Nexits)
        if verbose
            println("Solving logk $logk")
        end
        Nstart = Nexit - Nsubhorizon
        prob = ODEProblem(odefunc,Matrix{SolT}(I,2*d,2*d),(SolT(Nstart),SolT(Nend)),logk)
        psol = DifferentialEquations.solve(prob,options.scalars_solver,
        reltol=options.scalars_rtol,abstol=options.scalars_atol,
        maxiters=options.scalars_max_steps,save_everystep=options.save_everystep
        ,saveat=options.saveat)# automatically choose algorithm?
        G(N) = psol(N)
        push!(Gs,G)
    end
    if verbose
        println("Tracing out powerspectra")
    end
    Pi(N) = sol(N)[d+1:2*d]
    gNe(N) = vcat(-reshape(g(N),d,d)*Pi(N)/(2*Eh(N)),zeros((d)))
    estar = div(length(Nexits),2)+1
    Nstar = Nexits[estar]
    τ = -1/(exp(Nstar+logainit)*H(Nstar))
    Σ_BD = Array{SolT}(undef,length(logks),2*d,2*d)
    for (ek,logk) in enumerate(logks)
        Nstart = Nexits[ek] - Nsubhorizon
        kτ = -exp(logk - Nstart - logainit)/H(Nstart)
        gi = reshape(ginv(Nstart),d,d)
        S = 1/2 * H(Nstart)^2 * abs(kτ)^2 *
        [gi -gi;
         -gi gi.*abs(kτ)^2]
        Σ_BD[ek,:,:] = S
    end
    function Σ(ek,N)
        G = convert(Array{BigFloat},Gs[ek](N))
        S = convert(Array{BigFloat},Σ_BD[ek,:,:])
        @assert !any(isnan,S)
        @assert !any(isinf,S)
        @assert !any(isnan,G)
        @assert !any(isinf,G)
        sig = G*S*G'
        #sig = Array{Float64}(undef,2*d,2*d)
        #for a in 1:2*d
        #    for b in 1:2*d
        #        sig[a,b] = 0.0
        #        for c in 1:2*d
        #            for d in 1:2*d
        #                sig[a,b] += G[a,c]*G[b,d]*S[c,d]
        #                if isnan(sig[a,b])
        #                    println("$(G[a,c]), $(G[b,d]), $(S[c,d])")
        #                end
        #            end
        #        end
        #    end
        #end
        @assert !any(isnan,sig)
        return sig
    end
    function Pζ(ek,N)
        gN = convert(Array{BigFloat},gNe(N))
        @assert !any(isnan,gN)
        #if any([x<0 for x in gN])
        #    println(ek,N)
        #end
        #@assert(all([x>0 for x in gN]),"negative in gN")
        S = convert(Array{BigFloat},Σ(ek,N))
        @assert !any(isnan,S)
        #@assert(all([x>0 for x in S]),"negative in Σ")
        #if ek == estar
        #    println(Gs[ek](N))
        #    println(gN)
        #    println(S)
        #end
        P = gN'*S*gN
        if isnan(P)
            println(gN)
            println(S)
        end
        @assert !isnan(P)
        #P = 0.0
        #for a in 1:2*d
        #    for b in 1:2*d
        #        P += gN[a]*gN[b]*S[a,b]
        #    end
        #end
        P /= (2*pi^2)
        #@assert(P>0,"negative P")
        return P
    end

    # adapted for curved space -- TODO need a unit test
    function Piso(ek,N)
        gN = gNe(N)
        gi = reshape(ginv(N),d,d)
        viα = LinearAlgebra.nullspace(convert(Array{SolT,2},(gi*gN[1:d])'))
        S = Σ(ek,N)
        P = 0.0
        for a in 1:d
            for b in 1:d
                for i in 1:(d-1)
                    P += viα[a,i]*S[a,b]*viα[b,i]
                end
            end
        end
        P *= 1/(2*Eh(N))
        P /= (2*pi^2)
        #if verbose
            #@assert(P>0,"negative P")
        #end
        return P
    end

    Pζs = BigFloat[]
    #Pζs_hc = Float64[]
    Pisos = BigFloat[]
    for (ek,logk) in enumerate(logks)
        append!(Pζs,Pζ(ek,Nend))
        if d>1
            append!(Pisos,Piso(ek,Nend))
        end
        #append!(Pζs_hc,Pζ(ek,Nhc+3.0))
    end

    # fit for ns, as
    # logPζ = As + (ns-1) log(k/k⋆) + 1/2 αs log(k/k⋆)^2 + ...
    if verbose
        println("Fitting powerspectra")
    end
    p = poly_fit(logks.-logks[estar],log.(Pζs),2)
    As = exp(p[1])
    ns = 1+p[2]
    as = p[3]*2
    if verbose
        println("As = $As")
        println("ns = $ns")
        println("αs = $as")
    end
    #p = poly_fit(logks.-logks[estar],log.(Pζs_hc),2)
    #Ashc = exp(p[1])
    #nshc = 1+p[2]
    #ashc = p[3]*2
    #if verbose
    #    println("Ashc = $Ashc")
    #    println("nshc = $nshc")
    #    println("αshc = $ashc")
    #end
    #println(errors)
    if options.compute_isocurvature && d > 1
        p = poly_fit(logks.-logks[estar],log.(Pisos),2)
        Asiso = exp(p[1])
        nsiso = 1+p[2]
        asiso = p[3]*2
        if verbose
            println("Asiso = $Asiso")
            println("nsiso = $nsiso")
            println("αsiso = $asiso")
        end
    end

    Ns = (Nend-Nexit_to_end-Nsubhorizon):0.1:Nend
    rs = Array{BigFloat}(undef,length(Ns))
    if options.compute_tensors
        global rs
        if verbose
            println("Computing tensor perturbations at kstar")
        end
        global r
        r = NaN
        tensor_eom = funcs["Tensor_eom"]
        #tensor_jac = funcs["Tensor_jac"]

        #function Gamjac(J,Gam,logk,t)
        #    J[:,:] = reshape(tensor_jac(t,logainit,logk,Gam[1,1],Gam[1,2],Gam[2,2],[sol(t)[1:2*d];pvals]...),4,4)
        #end
        function Gamprime(dGam,Gam,logk,t)
            dGam[:,:] = reshape(tensor_eom(t,logainit,logk,Gam[1,1],Gam[1,2],Gam[2,2],[sol(t)[1:2*d];pvals]...),2,2)
        end
        Nexit = Nexits[estar]
        logk = logks[estar]
        if verbose
            println("Solving logk $logk")
        end
        Nstart = Nexit - Nsubhorizon
        kτH = -exp(logk - Nstart - logainit)
        kτ = kτH/H(Nstart)
        GamIC = [ 1 -1 ; -1  abs(kτ)^2] * abs(kτH)^2
        odefunc = ODEFunction(Gamprime)#,jac=Gamjac)
        tprob = ODEProblem(odefunc,GamIC,(Nstart,Nend),logk)
        tsol = DifferentialEquations.solve(tprob,options.tensors_solver,
        reltol=options.tensors_rtol,abstol=options.tensors_atol,
        maxiters=options.tensors_max_steps,save_everystep=options.save_everystep,
        saveat=options.saveat)# automatically choose algorithm?
        Gam(N) = tsol(N)
        #push!(Gs,G)
        r = 4*Gam(Nend)[1,1]/As / (2*pi^2)
        if verbose
            println("Tensor to scalar ratio: $r")
        end
        rf(N) = 4*Gam(N)[1,1]/As / (2*pi^2)
        rs = rf.(Ns)
    else
        global rs
        global r
        r = NaN
        rs = []
    end


    PζsN = [Pζ(estar,N) for N in Ns]
    PisosN = [Piso(estar,N) for N in Ns]

    tsol = Dict()
    tsol["As"] = As
    tsol["ns"] = ns
    tsol["as"] = as
    tsol["r"] = r
    tsol["PzN"]=PζsN
    tsol["PisoN"]=PisosN
    tsol["rN"] = rs
    tsol["Pzk"]=Pζs
    tsol["Pisok"]=Pisos
    tsol["lnk"] = logks.-logks[estar]
    tsol["pert_Ne"]=Ns

    return tsol
end
