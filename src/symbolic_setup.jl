function inflation_var_setup(d)
    f = [ symbols("ϕ$i",real=true) for i in 1:d ]
    fp= [ symbols("π$i",real=true) for i in 1:d ]
    return (f,fp)
end

function sum_separableQ(d,gv,f)
    for i in 1:d
        for j in 1:d
            if i != j && occursin(string(convert(Expr,f[i])),string(convert(Expr,gv[j])))
                return false
            end
        end
    end
    return true
end
function make_funcs(quotes)
    funcs = Dict{String,Function}()
    for (k,v) in quotes
        funcs[k] = @RuntimeGeneratedFunction(v)
    end
    return funcs
end

function inflation_setup(d,V,G,params=[];options::SetupOptions = SetupOptions())
    apply_hca = false
    f,fp = inflation_var_setup(d)
    v = V(f,params)
    g = G(f,params)
    if options.allow_caching
        hashed = hash((d,Quote(v),Quote.(g),Quote.(params),options.background_only,cache_version))
        fname = "funcs$hashed.bin"
        println("Looking for cached symbolic functions in $fname")
        if isfile(fname)
            println("Loading cached symbolic ϵV, Mab, etc...")
            open(fname,"r") do f
                quotes = deserialize(f)
            end
            if options.return_quotes
                return quotes
            else
                funcs = make_funcs(quotes)
                return funcs
            end
        end
    end
    println("Solving for symbolic ϵV, Mab, etc and compiling...")
    if d == 1
        ginv = [g[1,1]^(-1)]
    else
        ginv = options.Ginv(G,f,params)
    end
    if options.simplification > 1
        println("simplifying g, ginv")
        g = SymPy.simplify.(g)
        ginv = SymPy.simplify.(ginv)
    end
    gv = Sym[diff(v,f[i]) for i in 1:d ]
    if options.attempt_HCA
        apply_hca = sum_separableQ(d,gv,f)
    end
    if apply_hca
        println("Detected sum-separable potential, applying HCA")
    end
    ggv = [diff(v, f[i], f[j]) for i in 1:d, j in 1:d]
    if options.simplification >= 2
        println("simplifying ggv")
        ggv = simp.(ggv)
    end
    ev = epsv(gv,ginv,v)
    eh = epsh(fp,g)
    hub = hubble(v,eh)
    et = eta(fp,eh,hub,gv)
    om = omega(eh,ev,et)
    om3 = om
    if d > 1
        om3 = omega3(fp,gv,g,ginv,hub,eh)
    end
    println("starting h")
    #h,hl = options.christoffel(f,g,ginv)
    h,hl = options.christoffel(f,g,ginv)
    println("starting pi_eom")
    pi_eom = make_pi_eom(fp,h,ginv,gv,hub,eh;simplification=options.simplification)
    println(sympy.count_ops.(pi_eom))
    cse_pi_eom = sympy.cse(pi_eom)
    println(sympy.count_ops.(cse_pi_eom))
    if !options.background_only
        println("starting Rpp")
        Rpp = Rππ_fast(f,fp,h,hl,g,ginv,simplification=options.simplification)
        println(sympy.count_ops.(Rpp))
        #println(sympy.count_ops.(simplify(pi_eom[1])))
        geo_eom = make_geodesic_eom(fp,h)
        println("starting mass matrix")
        #pi_jac = make_pi_jac(f,fp,hub,eh,pi_eom)
        Covggv = covggv(ggv,h,gv,options.simp,options.simplification)
        println(sympy.count_ops.(Covggv))
        Mab = mab(Covggv,Rpp,fp,hub,gv,eh,g,ginv,options.simp,options.simplification)
        println(sympy.count_ops.(Mab))
        cse_mab = sympy.cse(Mab)
        println("cse Mab")
        println(sympy.count_ops.(cse_mab))

    #Sig = sig(fp,g)
    #S = s(ginv,g,Sig)
    #om2,S2 = omega2(fp,g,hub,h,pi_eom,Sig)
    #Mσσ,Mss = masses(g,Mab,Sig,S)
        println("done")
        println("making pert eom")

        loga = symbols("loga",real=true)
        logk = symbols("logk",real=true)
        t = symbols("t",positive=true)
        Gam = Array{Sym}(undef,2*d,2*d)
        for i in 1:2*d
            for j in 1:2*d
                Gam[i,j] = symbols("g$(i)o$j",real=true)
            end
        end
        GT = Array{Sym}(undef,2,2)
        for i in 1:2
            for j in 1:i
                GT[i,j] = symbols("gt$(i)o$j",real=true)
                GT[j,i] = GT[i,j]
            end
        end

        #pert_eom = make_pert_eom(t,Gam,fp,Mab,h,eh,loga,logk,hub)
        pert_uab = make_pert_uab(t,fp,Mab,h,eh,loga,logk,hub)
        pert_curv = make_pert_curv(fp,h)
        #pert_jac = make_pert_jac(Gam,fp,pert_eom)
        tensor_eom = make_tensor_eom(t,GT,eh,loga,logk,hub)

        #println(sympy.count_ops.(pert_eom))
        #println("cse pert_eom")
        #cse_pert_eom = sympy.cse(pert_eom)
        #println(sympy.count_ops.(cse_pert_eom[2][end]))
        #println("done")
        println("pert uab")
        println(sympy.count_ops.(pert_uab))
        println("cse pert_uab")
        cse_pert_uab = sympy.cse(pert_uab)
        println(sympy.count_ops.(cse_pert_uab[2][end]))
        println("done")
        println("cse pert_curv")
        cse_pert_curv = sympy.cse(pert_curv)
        println(sympy.count_ops.(cse_pert_curv[2][end]))
        println("done")
        println("tensor eom ops:")
        println(sympy.count_ops.(tensor_eom))

        println("done making pert eom")
        #tensor_jac = make_tensor_jac(GT,tensor_eom)

        #fd = hub*fp
        #if d ==1
        #    sig = sqrt((fd*g*fd)[1,1])
        #else
        #    sig = sqrt((fd'*g*fd)[1,1])
        #end
        #t = fd / sig
        #w = Array{Expr}(undef,d)
        #for a in 1:d
        #    w[a] = 0
        #    for b in 1:d
        #        w[a] += 1/sig*(ginv[a,b]-t[a]*t[b])*gv[b]
        #    end
        #end
        #if d == 1
        #    s = w / sqrt((w*g*w)[1,1])
        #else
        #    s = w / sqrt((w'*g*w)[1,1])
        #end
        #
    end
    if apply_hca
        va = [ integrate(gv[i],f[i]) for i in 1:d]
        va .+= (v - sum(va))/d
        ea = [gvi^2 / v^2 / 2 for gvi in gv ]
        etaa = [ggvi / v for ggvi in diag(ggv)]
        Pz_hca = pz_hca(v,ea,va)
        Ns_hca = ns_hca(v,ea,etaa,eh,va)
        println("Analytic Pζ (HCA):\n $Pz_hca")
        println("Analytic ns (HCA):\n $Ns_hca")
    end
    println("quoting")
    vars = [f;fp;params]
    quotes = Dict{String,Any}()
    quotes["Eh"] = QuoteFn("_Eh",eh,vars)
    println("pi_eom")
    if d > 1
        quotes["Pi_eom"] = QuoteFnArrCSE("_Pieom",cse_pi_eom,vars)
    else
        quotes["Pi_eom"] = QuoteFnArr("_Pieom",pi_eom,vars)
    end
    quotes["V"] = QuoteFn("_Pot",v,vars)
    println("ev")
    quotes["Ev"] = QuoteFn("_Ev",ev,vars)
    quotes["Et"] = QuoteFn("_Et",et,vars)
    println("om")
    quotes["Om"] = QuoteFn("_Om",om,vars)
    #sympy.count_ops(om3)
    cse_om3 = sympy.cse(om3)
    #sympy.count_ops.(cse_om3)
    quotes["OmAlt"] = QuoteFnCSE("_Om3",cse_om3,vars)
    quotes["H"] = QuoteFn("_H",hub,vars)
    println("g")
    quotes["ginv"] = QuoteFnArr("_ginv",Sym[ginv...],vars)
    quotes["g"] = QuoteFnArr("_g",Sym[g...],vars)
    if !options.background_only
        println("ggv")
        quotes["Vab"] = QuoteFnArr("_ggv",ggv,vars)
        quotes["Geo_eom"] = QuoteFnArr("_Geoeom",geo_eom,vars)
        println("mab")
        quotes["Mab"] = QuoteFnArrCSE("_Mab",cse_mab,vars)
        #mσσquote = QuoteFn(Mσσ,vars)
        #mssquote = QuoteFn(Mss,vars)
        #chquote = :(function Chr($(vars...)) return [$(h...)] end)
        println("pert eom")
        #peomquote = QuoteFnArr("_Perteom",[pert_eom...],[t,loga,logk,Gam...,vars...])
        #peomquote = QuoteFnArrCSE("_Perteom",cse_pert_eom,[t,loga,logk,Gam...,vars...])
        quotes["Pert_uab"] = QuoteFnArrCSE("_Pertuab",cse_pert_uab,[t,loga,logk,vars...])
        quotes["Pert_curv"] = QuoteFnArrCSE("_Pertcurv",cse_pert_curv,[vars...])
        #pjacquote = :(function Pertjac($t,$loga,$logk,$(Gam...),$(vars...)) return [$(pert_jac...)] end )
        quotes["Tensor_eom"] = QuoteFnArr("_Tensoreom",[tensor_eom...],[t,loga,logk,unique(GT)...,vars...])
    end
    if apply_hca
        quotes["Pz_hca"] = QuoteFn("_Pzhca",Pz_hca,vars)
        quotes["ns_hca"] = QuoteFn("_Nshca",Ns_hca,vars)
    end
    #quotes = [ehquote,eomquote,vquote,ggvquote,evquote,etquote,omquote,hquote,mabquote,mσσquote,mssquote,ginvquote,gquote,peomquote,teomquote,pzhcaquote,nshcaquote]
    println("done")
    if options.allow_caching
        open(fname,"w") do f
            serialize(f,quotes)
        end
    end
    if options.return_quotes
        return quotes
    else
        lambdified = make_funcs(quotes)
        return lambdified
    end
end


function make_pi_eom(Pi,h,ginv,va,H,e;simplification=0)
    d = length(Pi)
    curv = Array{Any}(undef,d)
    vupa = Array{Any}(undef,d)
    preeom = Array{Any}(undef,d)
    eom = Array{Any}(undef,d)
    for i in 1:d
        println(i)
        curv[i] = 0
        vupa[i] = 0
        preeom[i] = 0
        for j in 1:d
            vupa[i] += ginv[i,j]*va[j]/(H^2)
            for k in 1:d
                curv[i] += h[i,j,k]*(Pi[j]*Pi[k])
            end
        end
        preeom[i] = vupa[i]+curv[i]
        if simplification >= 3
            preeom[i] = simp(vupa[i]+curv[i])
        end
    end
    for i in 1:d
        eom[i] = -preeom[i]
        eom[i] -= (3-e)*Pi[i]
    end
    #println([-( (6-:p1^2)/:f1 + (3-:p1^2 / 2) * :p1 )])
    #return [-( (6-:p1^2)/:f1 + (3-:p1^2 / 2) * :p1 )]
    return eom
end

function make_geodesic_eom(Pi,h)
    d = length(Pi)
    eom = Array{Sym}(undef,d)
    for i in 1:d
        eom[i] = -Pi'*h[i,:,:]*Pi
    end
    eom
end

function make_pi_jac(Phi,Pi,H,eh,eom)
    vars = [Phi;Pi]#;H]
    d = length(Phi)
    l = length(vars)
    J = Array{Any}(undef,l,l)
    J[1:d,1:d] = zeros(d,d)
    J[1:d,d+1:2*d] = Matrix{Float64}(I,d,d)
    #J[end,1:d] = zeros(d)
    #J[end,d+1:2*d] = zeros(d)
    #J[end,end] = -eh
    for i in d+1:2*d
        for j in 1:l
            J[i,j] = df(eom[i-d],vars[j])
        end
    end
    return J
end

function make_pert_eom(t,G,p,mab,h,eh,loga,logk,hub)
    d = length(p)
    u = Array{Sym}(undef,2*d,2*d)
    curv = Matrix{Sym}(undef,2*d,2*d)
    u[1:d,1:d] = zeros((d,d))
    u[1:d,d+1:end] = Matrix{Float64}(I,d,d)
    #u[d+1:end,1:d] = -(Matrix{Int64}(I,d,d)*exp(2*(logk-t-loga-log(hub)))) - mab
    #u[d+1:end,d+1:end] = Matrix{Int64}(I,d,d)*(eh-3)
    for i in d+1:2*d
        for j in 1:d
            u[i,j] = -mab[i-d,j]
            if i-d == j
                u[i,j] -= exp(2*(logk-t-loga-log(hub)))
            end
        end
        for j in d+1:2*d
            u[i,j] = 0
            if i == j
                u[i,j] = eh-3
            end
        end
    end
    for i in 1:2*d, j in 1:2*d
        curv[i,j] = 0
        for m in 1:d, n in 1:d
            if i <= d
                curv[i,j] += h[i,m,n]*G[m,j]*p[n]
            else
                curv[i,j] += h[i-d,m,n]*G[m+d,j]*p[n]
            end
        #    #if j <= d
        #    #    curv[i,j] += h[j,m,n]*G[i,m]*p[n]
        #    #else
        #    #    curv[i,j] += h[j-d,m,n]*G[i,m+d]*p[n]
        #    #end
        end
        if simplification >= 2
            curv[i,j] = simp(curv[i,j])
            println(i,j)
        end
    end
    #println(u)
    #println(curv)
    return u*G .- curv
end

function make_pert_uab(t,p,mab,h,eh,loga,logk,hub)
    d = length(p)
    u = Array{Sym}(undef,2*d,2*d)
    curv = Matrix{Sym}(undef,2*d,2*d)
    u[1:d,1:d] = zeros((d,d))
    u[1:d,d+1:end] = Matrix{Float64}(I,d,d)
    u[d+1:end,1:d] = -(Matrix{Int64}(I,d,d)*exp(2*(logk-t-loga-log(hub)))) - mab
    u[d+1:end,d+1:end] = Matrix{Int64}(I,d,d)*(eh-3)

    return u
end
function make_pert_curv(p,h)
    d = length(p)
    hp = Array{Sym}(undef,d,d)
    for i in 1:d, j in 1:d
        hp[i,j] = sum([h[i,j,k]*p[k] for k in 1:d])
    end
    return hp
end

function make_pert_jac(G,p,pert_eom)
    d = length(p)
    vars = [G...]
    eom = [pert_eom...]
    l = length(vars)
    J = Array{Any}(undef,l,l)
    for i in 1:l
        for j in 1:l
            J[i,j] = df(eom[i],vars[j])
        end
    end
    return J
end

function make_tensor_eom(t,G,eh,loga,logk,hub)
    w = [ 0 1; -exp(2*(logk-t-loga-log(hub))) (eh-3)]
    dG = [ sum(w[a,c]*G[c,b] + w[b,c]*G[a,c] for c in 1:2) for a in 1:2 for b in 1:2]
    return dG
end

function make_tensor_jac(G,teom)
    vars = [G...]
    eom = [teom...]
    l = length(vars)
    J = Array{Any}(undef,l,l)
    for i in 1:l
        for j in 1:l
            J[i,j] = df(eom[i],vars[j])
        end
    end
    return J
end

function hubble(v,e)
    return sqrt(v/(3-e))
end
#function chrπ(Phi,g,ginv=inv(g))
#    """
#    returns the matrix ``Γ^a_{bc} π^c``, since this can be more efficiently calculated
#    than the christoffel symbols themselves
#    """
#    d = length(Pi)
#    gp = Array{Sym}(undef,d,d)
#    for j in 1:d, k in 1:d
#        gp[j,k] = 0
#        for l in 1:d, m in 1:d
#            gp[j,k] += 1/2*(ginv[j,m]*(diff(g[m,l],Phi[k])+diff(g[k,m],Phi[l])-diff(g[k,l],Phi[m])))*Pi[l]
#        end
#    end
#    return gpp
#
#end
#
#function chrππ(Phi,Pi,g,ginv=inv(g))
#    """
#    returns the vector ``Γ^a_{bc} π^b π^c``, since this can be more efficiently calculated
#    than the christoffel symbols themselves
#    """
#    d = length(Pi)
#    gpp = Array{Sym}(undef,d)
#    for j in 1:d
#        gpp[j] = 0
#        for k in 1:d, l in 1:d, m in 1:d
#            gpp[j] += 1/2*(ginv[j,m]*(diff(g[m,l],Phi[k])+diff(g[k,m],Phi[l])-diff(g[k,l],Phi[m])))*Pi[k]*Pi[l]
#        end
#    end
#    return gpp
#end

function christoffel(Phi,g,ginv=inv(g),simplification=1)
    """
    returns christoffel symbols as ``Γ^a_{bc}``
    """
    d = length(Phi)
    h  = Array{Sym}(undef,d,d,d)
    hl = Array{Sym}(undef,d,d,d)
    for j in 1:d, k in 1:d, l in k:d
        h[j,k,l]=0
        println(j,k,l)
        for m in 1:d
            h[j,k,l] += 1/2*(ginv[j,m]*(diff(g[m,l],Phi[k])+diff(g[k,m],Phi[l])-diff(g[k,l],Phi[m]) ))
        end
        hl[j,k,l] = 1/2*(diff(g[j,l],Phi[k])+diff(g[k,j],Phi[l])-diff(g[k,l],Phi[j]) )
        hl[j,k,l] = cancel(hl[j,k,l])
        h[j,k,l]  = cancel(h[j,k,l])
        h[j,l,k]  = h[j,k,l]
        hl[j,l,k] = hl[j,k,l]
    end
    return h,hl
end

function Rππ_fast(Phi,Pi,h,hl,g,ginv;simplification=1)
    """
    compute a common contraction of the Riemann curvature tensor ``R^α_{βγδ} π^β π^γ``
    try to use some of the riemann tensor symmetries to help, internally use R_{αβγδ}
    """
    function copy_symm!(R,i,j,k,l)
        # antisymm on first and last pair, and symm on swapping pairs
        R[i,j,l,k] = -R[i,j,k,l]
        R[j,i,k,l] = -R[i,j,k,l]
        R[j,i,l,k] = R[i,j,k,l]
        R[k,l,i,j] = R[i,j,k,l]
        R[l,k,i,j] = -R[i,j,k,l]
        R[k,l,j,i] = -R[i,j,k,l]
        R[l,k,j,i] = R[i,j,k,l]
    end

    d = length(Phi)
    R = Array{Sym}(undef,d,d,d,d)
    for i in 1:d
        R[i,i,:,:] = zeros(Int,d,d)
        R[:,:,i,i] = zeros(Int,d,d)
    end
    iterations = 0
    for i in 1:d, j in 1:d, k in 1:d, l in 1:d
        if !isassigned(R,i,j,k,l)
            # use lower index version for nice symmetries
            R[i,j,k,l] = -(diff(hl[l,j,k],Phi[i]) - diff(hl[l,i,k],Phi[j]))
            R[i,j,k,l] -= sum([(h[p,i,k]*hl[p,j,l] - h[p,j,k]*hl[p,i,l]) for p in 1:d])
            if simplification >= 2
                R[i,j,k,l] = simp(R[i,j,k,l])
            end
            # take care of (skew-)symmetries
            copy_symm!(R,i,j,k,l)
            #first Bianchi identity
            if !isassigned(R,k,l,i,j) && isassigned(R,j,k,l,i)
                R[k,l,i,j] = -R[i,j,k,l] - R[j,k,l,i]
                copy_symm!(R,k,l,i,j)
            elseif !isassigned(R,j,k,l,i) && isassigned(R,k,l,i,j)
                R[j,k,l,i] = -R[i,j,k,l] - R[k,l,i,j]
                copy_symm!(R,j,k,l,i)
            end
            iterations += 1
            println(iterations) # few over ideal, need to check correctness
        end
    end
    println("calculating Rππ")
    Rpp = Array{Sym}(undef,d,d)
    for i in 1:d, l in 1:d
        Rpp[i,l] = 0
        for j in 1:d, k in 1:d, m in 1:d
            Rpp[i,l] += ginv[i,m]*R[m,j,k,l]*Pi[j]*Pi[k]
        end
        if simplification >= 2
            Rpp[i,l] = simp(Rpp[i,l])
        end
    end
    return Rpp
end

function riemann(Phi,h)
    """
    inefficiently compute Riemann curvature tensor as ``R^α_{βγδ}``
    """
    d = length(Phi)
    R = Array{Any}(undef,d,d,d,d)
    for i in 1:d, j in 1:d, k in 1:d, l in 1:d
        R[i,j,k,l] = diff(h[i,l,j],Phi[k]) - diff(h[i,k,j],Phi[l])
        for m in 1:d
            R[i,j,k,l] += h[i,k,m]*h[m,l,j] - h[i,l,m]*h[m,k,j]
        end
    end
    return R
end

function ricci_tensor(R)
    """
    compute Ricci tensor ``R_{μν} = R^α_{μαν}``
    """
    Ric = Array{Expr}(0,d,d)
    q = CartesianIndices(Ric)
    for (i,j) in q
        for k in 1:d
            Ric[i,j] += R[k,i,k,j]
        end
    end
    return Ric
end
function ricci_scalar(Ric,ginv)
    """
    compute Ricci scalar ``R = g^{μν} R_{μν}``
    """
    return tr(ginv*Ric)
end

function epsh(Pi,g)
    e = 0
    d = length(Pi)
    for i in 1:d
        for j in 1:d
            e += Pi[i] * g[i,j] * Pi[j] / 2
        end
    end
    return e
end

function epsv(gv,ginv,v)
    d = length(gv)
    ev = 0
    for i in 1:d
        for j in 1:d
            ev += ginv[i,j]*(gv[i]*gv[j]) / (2*v^2)
        end
    end
    return ev
end

function eta(Pi,e,H,va)
    d = length(Pi)
    vp = 0
    for i in 1:d
        vp += va[i]*Pi[i]
    end
    return -6 + 2*e - vp/(H^2 * e)
end

function omega(e,ev,eta)
    return -e^2 - 6*ev + 9*ev/e - 1/4*(6 + eta)^2 + e*(6 + ev + eta)
end

function omega2(fp,g,H,h,pi_eom,sig)
    d = length(pi_eom)
    Σ(x) = reduce(Algebra.:+,x)
    P(x) = reduce(Algebra.:*,x)
    fpdot = pi_eom[:] # careful with the copy here, missing this is catastrophic
    signorm = 0
    for i in 1:d
        for j in 1:d
            signorm += sig[i] * g[i,j] * sig[j]
        end
    end
    signorm = sqrt(signorm)
    t = sig ./ 2
    w = fpdot
    w1 = t*(Σ([P([fp[i],fpdot[j],g[i,j]]) + P([fpdot[i],fp[j],g[i,j]]) for i in 1:d for j in 1:d]))
    w2 = t*Σ([P([fp[i],fp[j],fp[k],(h[l,k,i]*g[l,j]+h[l,k,j]*g[i,l])]) for i in 1:d for j in 1:d for k in 1:d for l in 1:d])
    w .-= w2
    w .-= w1
    w *= (H/signorm)
    Om = 0
    wnorm = 0
    for i in 1:d
        for j in 1:d
            wnorm += w[i]*g[i,j]*w[j]
        end
    end
    Om = wnorm / H^2
    wnorm = (wnorm)^(:q)
    s = [wa / wnorm for wa in w]
    s = [sub("q=1/2",sa) for sa in s]
    return Om,s
end

function omega3(fp,gv,g,ginv,hub,eh)
    vu = ginv*gv
    vu = vu ./ sqrt(vu'*g*vu)
    Vv = sqrt(gv'*ginv*gv)
    perpuu = ginv .- vu'.*vu
    piperp = perpuu*g*fp
    pnorm = sqrt(piperp'*g*piperp)
    return (pnorm * Vv / (2*eh*hub^2))^2
end
function covggv(ggv,h,gv,simp,simplification)
    d = length(gv)
    Covggv = Array{Sym}(undef,d,d)
    for i in 1:d
        for j in i:d
            Covggv[i,j] = ggv[i,j]
            for k in 1:d
                Covggv[i,j] -= h[k,i,j]*gv[k]
            end
            if simplification >= 1
                Covggv[i,j] = simp(Covggv[i,j])
            end
            Covggv[j,i] = Covggv[i,j]
        end
    end
    return Covggv
end

function mab(Covggv,Rpp,fp,hub,gv,eh,g,ginv,simp,simplification=1)
    d = length(gv)
    Mab = Array{Sym}(undef,d,d) # actually Mab/H^2
    sa = Array{Sym}(undef,d)
    for w in 1:d
        for b in 1:d
            for a in 1:d
                sa[a] = 0
                if ginv[w,a] != 0
                    s1 = ginv[w,a] * (Covggv[a,b]/hub^2)
                    s2 = -Rpp[w,b]
                    #s2 = -ginv[w,a]*sum([(g[a,m]*R[m,l,i,b])*(fp[l]*fp[i]) for l in 1:d for i in 1:d for m in 1:d])
                    s3 = -ginv[w,a]*(3+eh)*sum([g[a,m]*fp[m] for m in 1:d])*sum([g[b,n]*fp[n] for n in 1:d])
                    sum1 = sum([g[b,n]*fp[n] for n in 1:d])
                    sum2 = sum([g[a,m]*fp[m] for m in 1:d])
                    s4 = prod([ginv[w,a],2,eh,sum1,sum2])
                    s5 = -prod([ginv[w,a],sum2,(-gv[b]/hub^2 + sum1*(eh-3))])
                    s6 = -prod([ginv[w,a],sum1,(-gv[a]/hub^2 + sum2*(eh-3))])
                    if simplification >= 3
                        sa[a] = simp(s6+s5+s4+s3+s2+s1)
                    else
                        sa[a] = s6+s5+s4+s3+s2+s1
                    end
                end
            end
            Mab[w,b] = sum(sa)
        end
    end
    return Mab
end

function pz_hca(v,ea,va)
    # see 1303.3611 for these
    ua = [vai / v for vai in va]
    return v / (24 * pi^2) * sum((ua.^2 ./ ea))
end
function ns_hca(v,ea,etaa,eh,va)
    d = length(va)
    # see 1303.3611 for these
    ua = [vai / v for vai in va]
    num = (1-sum([etaa[i] * ua[i]^2 / (2*ea[i]) for i in 1:d]))
    den = sum([ua[i]^2 / ea[i] for i in 1:d])
    ns = 1 - 2*eh - 4*num/den
    return ns
end
function sig(fp,g)
    d = length(fp)
    norm = 0
    for i in 1:d
        for j in 1:d
            norm += fp[i]*g[i,j]*fp[j]
        end
    end
    norm = sqrt(norm)
    return [p / norm for p in fp]
end
# this is wrong, only works for 2-field
function s(ginv,g,sig)
    d = length(sig)
    sij = Array{Any}(undef,d,d)
    for i in 1:d
        for j in 1:d
            sij[i,j] = ginv[i,j] - sig[i]*sig[j]
        end
    end
    s = sij[:,d]
    norm = 0
    for i in 1:d
        for j in 1:d
            norm += s[i]*g[i,j]*s[j]
        end
    end
    return [si / norm for si in s]
end

function masses(g,Mab,Sig,S)
    d = length(Sig)
    Mσσ = 0
    for i in 1:d
        for j in 1:d
            for k in 1:d
                Mσσ += g[i,k]*Mab[k,j]*Sig[i]*Sig[j]
            end
        end
    end
    Mss = 0
    if d > 1
        for i in 1:d
            for j in 1:d
                for k in d:1
                    Mss += (g[i,k]*Mab[k,j])*S2[i]*S2[j]
                end
            end
        end
    end
    return Mσσ, Mss
end
