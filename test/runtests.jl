using Inflation
using SymPy
using LinearAlgebra
using Test

@testset "inflation models" begin
    @testset "single-field quadratic" begin
        d = 1
        function G(Phi,p)
            return [1.0]
        end

        m = symbols("m")
        params = [m]
        pvalues = [1.0e-10]
        function V(Phi,p)
            g = G(Phi,p)
            m = p[1]
            v = 0.0
            for i in 1:d
                for j in 1:d
                    v += (Phi[i] * m * Phi[j])/2.0
                end
            end
            return v
        end

        opt = SetupOptions(
        allow_caching = false
        )
        funcs = inflation_setup(d,V,G,params,options=opt)
        Phi0 = [18.0]
        Pi0  = [0.0]

        sol1 = background_evolve(Phi0,Pi0,pvalues,funcs,verbose=true)
        @test sol1.retcode == :Terminated
        @test sol1.t[end] ≈ 81.9 rtol=0.02
        eh = funcs["Eh"]([sol1.u[end];pvalues]...)
        @test eh ≈ 1.0 atol=1e-8
        tsol = transport_perturbations(sol1,pvalues,funcs,verbose=true)
        ns = tsol["ns"]
        as = tsol["as"]
        As = tsol["As"]
        r = tsol["r"]
        # results from the same model in Multi-Mode-Code
        MMC_Ne = 81.607534
        #MMC_al = 14.67
        MMC_As = 5.098295961e-9
        MMC_Aiso = 1.591443943
        MMC_Pent = 9.881312917
        MMC_Pcross = -10
        MMC_r = 0.1440518
        MMC_ns = 0.9633723789
        MMC_alpha_s = -6.854134226e-04
        MMC_f_NL = 2.842
        @test ns ≈ MMC_ns rtol=0.01
        @test as ≈ MMC_alpha_s rtol=0.05
        @test As ≈ MMC_As rtol=0.05
        @test r ≈ MMC_r rtol=0.02
    end
    @testset "single-field α-attractors" begin
        d = 1
        α = symbols("a",real=true)
        m = symbols("m",real=true)
        params = [α,m]
        pvalues = [1/2,1e-5]

        function G(Phi,p)
            return [6*α/(1-Phi[1]^2)^2]
        end
        function V(Phi,p)
            return α/2.0*(m^2*Phi[1]^2)
        end
        Phi0 = [.999]
        Pi0 = [0.0]

        opt = SetupOptions(
        allow_caching = false,
        attempt_HCA = false
        )

        funcs = inflation_setup(d,V,G,params,options=opt)
        sol1 = background_evolve(Phi0,Pi0,pvalues,funcs,verbose=true)

        #analytics
        Nexit_to_end = 55
        anal_ns = 1 - 2/Nexit_to_end
        anal_r  = 12*pvalues[1]/Nexit_to_end^2
        @test sol1.retcode == :Terminated
        @test sol1.t[end] ≈ 376.666 rtol=0.002
        tsol = transport_perturbations(sol1,pvalues,funcs,verbose=true)
        ns = tsol["ns"]
        as = tsol["as"]
        As = tsol["As"]
        r = tsol["r"]
        @test ns ≈ anal_ns rtol=0.02
        @test r ≈ anal_r rtol=0.08
        #@test as ≈ MMC_alpha_s rtol=0.05
        #@test As ≈ MMC_As rtol=0.05
    end
    @testset "double-field α-attractors" begin
        "see fig. 1 of 1803.09841 for this model"
        d = 2
        α,m,R = symbols("a m r",real=true)
        params = [α,m,R]
        pvalues = [1/6,1e-5,10]
        #α = 1/600
        function G(Phi,p)
            α = p[1]
            g = [ 6*α/(1-Phi[1]^2 - Phi[2]^2)^2 0; 0 6*α/(1-Phi[1]^2 - Phi[2]^2)^2]
            return g
        end
        #function Ginv(Phi,p)
        #    α = p[1]
        #    gi = [ :(1/(6*$α)*(1-$(Phi[1])^2 - $(Phi[2])^2)^2) :(0); :(0) :(1/(6*$α)*(1-$(Phi[1])^2 - $(Phi[2])^2)^2)]
        #    return gi
        #end
        function V(Phi,p)
            m = [p[2]^2 0.0; 0.0 p[3]*p[2]^2]
            v = 0
            for i in 1:d
                for j in 1:d
                    v += α/2.0*(m[i,j]*Phi[i]*Phi[j])
                end
            end
            return v
        end
        r0 = 0.999
        #r0 = 0.99
        θ0 = pi/4
        Phi0 = [r0*cos(θ0), r0*sin(θ0)]
        Pi0 = [0.0, 0.0]

        opt = SetupOptions(
        allow_caching = false,
        attempt_HCA = false
        )
        funcs = inflation_setup(d,V,G,params,options=opt)#,Ginv=Ginv)
        sol1 = background_evolve(Phi0,Pi0,pvalues,funcs,verbose=true)

        #analytics
        Nexit_to_end = 55
        Nang = 0
        Neff = Nexit_to_end - Nang
        anal_ns = 1 - 2/Neff
        anal_r  = 12*pvalues[1]/Neff^2
        @test sol1.retcode == :Terminated
        @test sol1.t[end] ≈ 129.567 rtol=0.002
        tsol = transport_perturbations(sol1,pvalues,funcs,verbose=true)
        ns = tsol["ns"]
        as = tsol["as"]
        As = tsol["As"]
        r = tsol["r"]
        @test ns ≈ anal_ns rtol=0.002
        @test r ≈ anal_r rtol=0.12 # this is about 10% higher than analytic estimate
        #@test as ≈ MMC_alpha_s rtol=0.05
        #@test As ≈ MMC_As rtol=0.05
    end
    @testset "2-field hyperinflation" begin
        d = 2
        L = symbols("l")
        params = [L]
        pvalues = [0.01]
        function G(Phi,p)
            #return [1 0;0 L^2*sinh(Phi[1]/L)^2]
            L = p[1]
            sh = sinh((Phi[1])/L)
            return [1 0;0 L^2*sh^2]
        end

        function V(Phi,p)
            L = p[1]
            return L^2 * Phi[1]^2 / 2
        end
        opt = SetupOptions(
        allow_caching = false,
        attempt_HCA = false
        )
        funcs = inflation_setup(d,V,G,params,options=opt)
        Phi0 = [0.6925,1e-24]
        Pi0 = [-0.0297881,4.8921e-29]

        sol2 = background_evolve(Phi0,Pi0,pvalues,funcs,verbose=true)
        @test sol2.retcode == :Terminated
        @test sol2.t[end] ≈ 23.05 rtol=0.02
        eh = funcs["Eh"]([sol2.u[end][1:end];pvalues]...)
        @test eh ≈ 1.0 atol=1e-8
    end
    @testset "2-field quadratic" begin
        d = 2
        params = [symbols(x,real=true,positive=true) for x in ["m","r"]]
        pvalues = [1e-5,10]
        function G(Phi,p)
            return Matrix{Int64}(I,d,d)
        end

        function V(Phi,p)
            m = [p[1]^2 0 ; 0 p[2]*p[1]^2]
            v = Sym(0.0)
            for i in 1:d
                for j in 1:d
                    v += (Phi[i] * m[i,j] * Phi[j])/2.0
                end
            end
            return v
        end

        Phi0 = [15.0,10.0]
        Pi0  = [0.0,0.0]
        opt = SetupOptions(
        allow_caching = false
        )
        funcs = inflation_setup(d,V,G,params,options=opt)
        sol2 = background_evolve(Phi0,Pi0,pvalues,funcs,verbose=true)
        @test sol2.retcode == :Terminated
        @test sol2.t[end] ≈ 82.407 rtol=0.02
        eh = funcs["Eh"]([sol2.u[end][1:end];pvalues]...)
        @test eh ≈ 1.0 atol=1e-8
        mT_Ne = 82.407
        mT_As = 9.8619e-9
        mT_ns = 0.919087
        mT_r = 0.144267
        tsol = transport_perturbations(sol2,pvalues,funcs,verbose=true)
        ns = tsol["ns"]
        as = tsol["as"]
        As = tsol["As"]
        r = tsol["r"]
        @test ns ≈ mT_ns rtol=0.02
        #@test as ≈ MMC_alpha_s rtol=0.05
        #MMC_alpha_s = -6.535329799e-04
        #MMC_f_NL = -1.5301e-02
        @test As ≈ mT_As rtol=0.05
        @test mT_r ≈ r rtol=0.02
    end
    @testset "5-field quadratic" begin
        d = 5
        params = [symbols("m$i",real=true,positive=true) for i in 1:d]
        pvalues = [1e-5,1e-5,sqrt(1e-9),1e-5,1e-4]
        function G(Phi,p)
            return Matrix{Int64}(I,d,d)
        end

        function V(Phi,p)
            m = p.^2
            v = 0.0
            for i in 1:d
                v += (Phi[i]^2 * m[i])/2.0
            end
            return v
        end

        Phi0 = [15.0,10.0,-5.0,10.0,3.0]
        Pi0  = [0.0,0.0,0.0,0.0,0.0]

        opt = SetupOptions(
        allow_caching = false
        )
        funcs = inflation_setup(d,V,G,params,options=opt)
        sol5 = background_evolve(Phi0,Pi0,pvalues,funcs,verbose=true)
        @test sol5.retcode == :Terminated
        @test sol5.t[end] ≈ 115.62897 rtol=0.02
        eh = funcs["Eh"]([sol5.u[end][1:end];pvalues]...)
        @test eh ≈ 1.0 atol=1e-8
        MMC_Ne = 1.1562897e+02
        #MMC_al = 14.67
        MMC_As = 5.105241913e-09
        MMC_Aiso = 4.217972911e-13
        MMC_Pent = 0
        MMC_Pcross = -2.026595982e-20
        MMC_r = 0.1440332439
        MMC_ns = 0.9631782590
        MMC_alpha_s = -6.535329799e-04
        MMC_f_NL = -1.5301e-02
        tsol = transport_perturbations(sol5,pvalues,funcs,verbose=true)
        ns = tsol["ns"]
        as = tsol["as"]
        As = tsol["As"]
        r = tsol["r"]
        @test ns ≈ MMC_ns rtol=0.01
        @test as ≈ MMC_alpha_s rtol=0.05
        @test As ≈ MMC_As rtol=0.05
        @test r ≈ MMC_r rtol=0.02
    end
end
