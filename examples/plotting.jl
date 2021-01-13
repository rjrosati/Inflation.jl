### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ b562e4c8-d339-11ea-3903-a936ef2608e6
using JLD2

# ╔═╡ e58fb87e-d339-11ea-15e9-755175e9db02
using Plots

# ╔═╡ 5f25df2c-d33e-11ea-2349-b3cde1716f8f
using LaTeXStrings

# ╔═╡ 9ba366a6-d33f-11ea-294e-57f4c58b40c0
md"# Load saved data"

# ╔═╡ e93aba66-d339-11ea-0296-b360f581dba3
begin
	data = Dict{String,Any}()
	jldopen( "Nflation_example.jld2", "r") do f
		for k in keys(f)
			data[k] = read(f,k)
		end
	end
	data
end

# ╔═╡ adca77d4-d33f-11ea-2a01-31ef55bf0781
md"# Background evolution"

# ╔═╡ 477e130e-d33b-11ea-2c9b-8300969f6553
begin
	l = @layout [a b ; c d]
	ehp = plot(data["back_Ne"],data["eh"], ylabel=L"\epsilon_H", xlabel=L"N_e")
	evp = plot(data["back_Ne"],data["ev"], ylabel=L"\epsilon_V", xlabel=L"N_e")
	etp = plot(data["back_Ne"],data["et"], ylabel=L"\eta_H", xlabel=L"N_e")
	omp = plot(data["back_Ne"],data["omega"], ylabel=L"\omega^2/H^2", xlabel=L"N_e")
	plot(ehp,evp,etp,omp,layout=l, legend=false, linewidth=2)
	#plot(ehp,evp,omp,layout=l, legend=false, linewidth=2)
end

# ╔═╡ a6dc05f4-d340-11ea-22a7-39b8356c7cd8
md"# Planck compatibility?"

# ╔═╡ 57e277c6-d345-11ea-2cc1-4fb4bfc0be6f
md"This is a rough plot based on the Planck 2018 data release, real planck data should be used if you want to be more careful."

# ╔═╡ b6b376d8-d340-11ea-2d9f-2f68ab41a13b
begin
	planck_ns = 0.9668
	planck_ns_std = 0.0037
	planck_r = 0
	planck_r_std = 0.058/1.95996
	
	num_sigma = 5
	function ellipseShape(h,k,a,b)
		θ = range(0,2π,length=300)
		h .+ a*cos.(θ), k .+ b*sin.(θ)
	end
	plot(ellipseShape(planck_ns,planck_r,planck_ns_std,planck_r_std), seriestype = [:shape,], lw=0.5, c = :red, linecolor=:black, legend=false, fillalpha = 0.2,label="Planck 2018 (up to $(num_sigma)σ)")
	for i in 2:num_sigma
		plot!(ellipseShape(planck_ns,planck_r,i*planck_ns_std,i*planck_r_std), seriestype = [:shape,], lw=0.5, c = :red, linecolor=:black, legend=false, fillalpha = 0.2,label="")
	end
	scatter!([data["ns"]],[data["r"]],color=:green,leg=true,label="This model")
	plot!(ylim=[0,:auto])
end

# ╔═╡ a2f15f0e-d352-11ea-2f7b-e16f13cbdd76
md"# Sub and super-horizon evolution of 2-pt functions"

# ╔═╡ 79dabff8-d356-11ea-0363-d1d94cdbc5c1
md"The powerspectra are plotted starting when Bunch-Davies initial conditions are imposed, and horizon exit is marked as a vertical red line"

# ╔═╡ b6bb0206-d352-11ea-3110-63c5067874a8
begin
	Nend = data["pert_Ne"][end]
	plot(data["pert_Ne"] .- Nend, data["PzN"],label=L"P_\zeta",lw=2)
	plot!(data["pert_Ne"] .- Nend, data["PisoN"],label=L"P_\mathrm{iso}",lw=2)
	plot!(data["pert_Ne"] .- Nend, data["rN"].*data["PzN"],label=L"P_\mathrm{tens}",lw=2)
	plot!(yscale=:log10)
	vline!([-55.0],color=:red,lw=2,linestyle=:dash,label="")
	plot!(xlabel=L"N_e",ylabel=L"P_{\circ}(k_\star,N)")
end

# ╔═╡ Cell order:
# ╠═b562e4c8-d339-11ea-3903-a936ef2608e6
# ╠═e58fb87e-d339-11ea-15e9-755175e9db02
# ╠═5f25df2c-d33e-11ea-2349-b3cde1716f8f
# ╟─9ba366a6-d33f-11ea-294e-57f4c58b40c0
# ╠═e93aba66-d339-11ea-0296-b360f581dba3
# ╟─adca77d4-d33f-11ea-2a01-31ef55bf0781
# ╠═477e130e-d33b-11ea-2c9b-8300969f6553
# ╟─a6dc05f4-d340-11ea-22a7-39b8356c7cd8
# ╟─57e277c6-d345-11ea-2cc1-4fb4bfc0be6f
# ╠═b6b376d8-d340-11ea-2d9f-2f68ab41a13b
# ╟─a2f15f0e-d352-11ea-2f7b-e16f13cbdd76
# ╟─79dabff8-d356-11ea-0363-d1d94cdbc5c1
# ╠═b6bb0206-d352-11ea-3110-63c5067874a8
