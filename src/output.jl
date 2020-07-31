function output_data(name,funcs,pvals,sol,tsol)
    Nend = sol.t[end]
    bNs = 0.0:0.01:Nend
    d = length(sol.u[1])÷2
    #Mσσ(N) = funcs["Mσσ"]([sol(N);pvals]...)
    #Mss(N) = funcs["Mss"]([sol(N);pvals]...)
    Eh(N) = funcs["Eh"]([sol(N)[1:2*d];pvals]...)
    Ev(N) = funcs["Ev"](sol(N)...,pvals...)
    Om(N) = funcs["Om"]([sol(N)[1:2*d];pvals]...)
    if haskey(funcs,"Mab")
        Mab(N) = funcs["Mab"]([sol(N)[1:2*d];pvals]...)
        mabs = Mab.(bNs)
    end
    #Mσσs = Mσσ.(bNs)
    #Msss = Mss.(bNs)
    ehs = Eh.(bNs)
    evs = Ev.(bNs)
    oms = Om.(bNs)
    phis = [u[1:d] for u in sol.(bNs)]
    pis = [u[d+1:end] for u in sol.(bNs)]
    jldopen(name,"w") do f
        write(f,"eh",ehs)
        write(f,"ev",evs)
        write(f,"omega",oms)
        if haskey(funcs,"Mab")
            write(f,"Mab",mabs)
        end
        write(f,"Phi",phis)
        write(f,"Pi",pis)
        write(f,"pvalues",pvals)
        write(f,"back_Ne",bNs)
        if !isnothing(tsol)
            for (key,val) in tsol
                write(f,key,val)
            end
        end
    end
end

function generate_default_plots(output_fname)
    # write to a pluto notebook
end
