using DelimitedFiles
using Plots

data = readdlm("results/temperature.txt")

r = data[:,1]
θ = data[:,2]
z = data[:,3]
T = data[:,4]

x = r .* cos.(θ)
y = r .* sin.(θ)

scatter(
    x, y, z,
    marker_z=T,
    markersize=2,
    color=:thermal,
    xlabel="x",
    ylabel="y",
    zlabel="z",
    title="Temperature Field T(r, θ, z)",
    legend=false
)
