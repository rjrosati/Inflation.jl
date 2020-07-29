@with_kw struct SetupOptions
    allow_caching::Bool = true
    return_quotes::Bool = false
    background_only::Bool = false
    Ginv::Function = (G,f,params) -> inv(G(f,params))
    christoffel::Function = christoffel
    simplification::Int = 1 # level of simplification, bigger is more aggressive
    # the simplification function, try simplify(x) for small Nf
    # factor(expand(x)) is more efficient at higher Nf
    # see https://docs.sympy.org/latest/tutorial/simplification.html for more
    simp::Function = x-> cancel(x)
    attempt_HCA::Bool = true # attempt to apply the horizon-crossing approximation
end
@with_kw struct BackgroundOptions
    background_atol::Float64 = 1e-14
    background_rtol::Float64 = 1e-14
    max_efolds::Float64 = 1e4
    min_efolds::Float64 = 0.0 # force eom evolution even if ϵ > 1 for this long
    max_steps::Float64 = Int(1e8)
    background_solver = Tsit5()
    end_epsilon::Float64 = 1.0 # the value of ϵ that ends inflation
    save_everystep::Bool = false
    saveat::Float64 = 1e-3
    Nexit_to_end::Float64 = 55
end
@with_kw struct PerturbationOptions
    # perturbations
    force_transport::Bool = true # compute perturbations even if inflation doesn't end
    compute_isocurvature::Bool = true # when d > 1, compute scalar powerspectra from all non-adiabatic 2-point functions
    compute_tensors::Bool = true # compute tensor perturbations
    Nsubhorizon::Float64 = 8. # number of efolds before horizon exit to impose Bunch-Davies
    Nkrange::Float64 = 5.0 # number of e-folds spanned by k-values
    Nexit_to_end::Float64 = 55
    kpivot::Float64 = 0.002 # Mpc^(-1)
    desired_ksteps::Int = 6 # odd number recommended, to measure ns right at pivot scale
    # scalar eom
    scalars_max_steps::Int = 4e8
    scalars_atol::Float64 = 1e-13
    scalars_rtol::Float64 = 1e-13
    # picking the right solver can be very problem dependent. https://docs.juliadiffeq.org/stable/tutorials/advanced_ode_example/#Choosing-a-Good-Solver-1
    # if Tsit5 fails, try lowering the atol and rtol and switching to Rodas5
    # All worth a try: Rosenbrock23, Rodas5, lsoda, TRBDF2, KenCarp4
    scalars_solver = Tsit5()
    save_everystep::Bool = false
    saveat::Float64 = 1e-1
    # tensor eom
    tensors_max_steps::Int = scalars_max_steps
    tensors_atol::Float64 = scalars_atol
    tensors_rtol::Float64 = scalars_rtol
    tensors_solver = scalars_solver
end
