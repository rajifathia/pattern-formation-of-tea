### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 22e611a8-9fe3-11f1-80f1-d533ea9fd220
begin
    using Pkg
    Pkg.activate("/root/pattern-formation-of-tea")

    using PlutoUI
    using WGLMakie
    WGLMakie.activate!()
    using DelimitedFiles
end

# ╔═╡ a70707e6-4fbf-40df-a08d-3dc88f30553b
Pkg.add("CairoMakie")

# ╔═╡ 10ef579a-daba-4ab6-a837-9a32adf6eaff
data = readdlm("results/temperature.txt")

# ╔═╡ d5c68057-3dbd-48d5-a88b-4dd8779a633e
begin
	r = data[:, 1]
	θ = data[:, 2]
	z = data[:, 3]
	T = data[:, 4]
end

# ╔═╡ 69b1b35e-3bed-4b1e-9e47-09c459fd645a
begin
    using CairoMakie
    
    fig = Figure(size = (800, 500))
    ax = Axis(
        fig[1, 1],
        xlabel = "z",
        ylabel = "Temperature",
        title = "Temperature Distribution"
    )
    
    lines!(ax, z, T)
    
    fig
end

# ╔═╡ b0c01d73-9c21-4827-83a6-61ff35b933d5
using Pkg

# ╔═╡ Cell order:
# ╠═22e611a8-9fe3-11f1-80f1-d533ea9fd220
# ╠═10ef579a-daba-4ab6-a837-9a32adf6eaff
# ╠═d5c68057-3dbd-48d5-a88b-4dd8779a633e
# ╠═a70707e6-4fbf-40df-a08d-3dc88f30553b
# ╠═b0c01d73-9c21-4827-83a6-61ff35b933d5
# ╠═69b1b35e-3bed-4b1e-9e47-09c459fd645a
