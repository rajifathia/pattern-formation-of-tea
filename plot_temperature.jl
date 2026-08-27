using DelimitedFiles
using Plots

data = readdlm("results/temperature.txt")

r = data[:,1]
theta = data[:,2]
z = data[:,3]
T = data[:,4]

# Take the θ = 0 cylindrical cross-section
mask = abs.(theta) .< 1e-10

r0 = r[mask]
z0 = z[mask]
T0 = T[mask]

# Unique grid coordinates
rs = sort(unique(r0))
zs = sort(unique(z0))

# Build temperature matrix T(r,z)
M = fill(NaN, length(zs), length(rs))

for n in eachindex(T0)
    i = findfirst(==(r0[n]), rs)
    k = findfirst(==(z0[n]), zs)
    M[k,i] = T0[n]
end

heatmap(
    rs,
    zs,
    M,
    xlabel = "Radial position r",
    ylabel = "Height z",
    title = "Temperature Field — θ = 0",
    color = :thermal,
    clims = (20,30),
    colorbar_title = "Temperature (°C)",
    aspect_ratio = :equal
)

savefig("results/temperature_rz_heatmap.png")

println("DONE: results/temperature_rz_heatmap.png")
println("Temperature range: ", minimum(T0), " to ", maximum(T0), " °C")
