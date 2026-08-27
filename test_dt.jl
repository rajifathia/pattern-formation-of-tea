\
 # ==========================================================
# BLOCK 1
# CYLINDRICAL GRID AND PARAMETERS
# ==========================================================

using LinearAlgebra

println("Initializing cylindrical computational domain...")

# ----------------------------------------------------------
# Geometry
# ----------------------------------------------------------

R = 1.0                 # Cylinder radius
H = 1.0                 # Cylinder height

# ----------------------------------------------------------
# Number of grid points
# ----------------------------------------------------------

Nr_points = 51          # Radial direction
Nθ_points = 64          # Azimuthal direction
Nz_points = 51          # Vertical direction

# ----------------------------------------------------------
# Grid spacing
# ----------------------------------------------------------

Δr = R/(Nr_points-1)
Δθ = 2π/Nθ_points
Δz = H/(Nz_points-1)

# ----------------------------------------------------------
# Coordinate arrays
# ----------------------------------------------------------

r = collect(0:Δr:R)

θ = collect(0:Δθ:(2π-Δθ))

z = collect(0:Δz:H)

# ----------------------------------------------------------
# Physical parameters
# ----------------------------------------------------------

ρ₀ = 1000.0             # Density (kg/m³)

ν = 1.0e-6              # Kinematic viscosity

κ = 1.4e-7              # Thermal diffusivity

D = 1.0e-9              # Concentration diffusivity

g = 9.81                # Gravity

β = 2.1e-4              # Thermal expansion coefficient

T₀ = 25.0               # Reference temperature

# ----------------------------------------------------------
# Time-stepping
# ----------------------------------------------------------

dt = 1.0e-3

Nt = 50

# ----------------------------------------------------------
# Boundary temperatures
# ----------------------------------------------------------

Tb = 30.0               # Bottom plate

Tt = 20.0               # Top plate

# ----------------------------------------------------------
# Allocate solution arrays
# ----------------------------------------------------------

ur = zeros(Nr_points,Nθ_points,Nz_points)

uθ = zeros(Nr_points,Nθ_points,Nz_points)

uz = zeros(Nr_points,Nθ_points,Nz_points)

p = zeros(Nr_points,Nθ_points,Nz_points)

T = fill(T₀,Nr_points,Nθ_points,Nz_points)

C = zeros(Nr_points,Nθ_points,Nz_points)

ur_star = zeros(Nr_points,Nθ_points,Nz_points)
uθ_star = zeros(Nr_points,Nθ_points,Nz_points)
uz_star = zeros(Nr_points,Nθ_points,Nz_points)

println("Cylindrical grid initialized successfully.")




 # ==========================================================
# BLOCK 2
# FINITE DIFFERENCE OPERATORS
# CYLINDRICAL COORDINATES
# ==========================================================

println("Loading cylindrical finite-difference operators...")

# ----------------------------------------------------------
# FIRST DERIVATIVES
# ----------------------------------------------------------

# ∂f/∂r
function dr(f,i,j,k)
    return (f[i+1,j,k] - f[i-1,j,k])/(2Δr)
end

# ∂f/∂θ (Periodic boundary)
function dθ(f,i,j,k)

    jp = (j == Nθ_points) ? 1 : j + 1
    jm = (j == 1) ? Nθ_points : j - 1

    return (f[i,jp,k] - f[i,jm,k])/(2Δθ)

end

# ∂f/∂z
function dz(f,i,j,k)
    return (f[i,j,k+1] - f[i,j,k-1])/(2Δz)
end

# ----------------------------------------------------------
# SECOND DERIVATIVES
# ----------------------------------------------------------

# ∂²f/∂r²
function drr(f,i,j,k)
    return (f[i+1,j,k] - 2f[i,j,k] + f[i-1,j,k])/(Δr^2)
end

# ∂²f/∂θ²
function dθθ(f,i,j,k)

    jp = (j == Nθ_points) ? 1 : j + 1
    jm = (j == 1) ? Nθ_points : j - 1

    return (f[i,jp,k] - 2f[i,j,k] + f[i,jm,k])/(Δθ^2)

end

# ∂²f/∂z²
function dzz(f,i,j,k)
    return (f[i,j,k+1] - 2f[i,j,k] + f[i,j,k-1])/(Δz^2)
end

# ----------------------------------------------------------
# CYLINDRICAL LAPLACIAN
#
# ∇²f =
# ∂²f/∂r²
# + (1/r)∂f/∂r
# + (1/r²)∂²f/∂θ²
# + ∂²f/∂z²
# ----------------------------------------------------------

function laplacian(f,i,j,k)

    ri = max(r[i], Δr/2)

    return drr(f,i,j,k) +
           (1/ri)*dr(f,i,j,k) +
           (1/ri^2)*dθθ(f,i,j,k) +
           dzz(f,i,j,k)

end

# ----------------------------------------------------------
# CYLINDRICAL ADVECTION
#
# u·∇f
# = ur ∂f/∂r
# + (uθ/r) ∂f/∂θ
# + uz ∂f/∂z
# ----------------------------------------------------------

function advection(f,ur,uθ,uz,i,j,k)

    ri = max(r[i], Δr/2)

    return ur[i,j,k]*dr(f,i,j,k) +
           (uθ[i,j,k]/ri)*dθ(f,i,j,k) +
           uz[i,j,k]*dz(f,i,j,k)

end

# ----------------------------------------------------------
# CYLINDRICAL DIVERGENCE
#
# ∇·u =
# (1/r)∂(rur)/∂r
# + (1/r)∂uθ/∂θ
# + ∂uz/∂z
# ----------------------------------------------------------

function divergence(ur,uθ,uz,i,j,k)

    ri = max(r[i], Δr/2)

    jp = (j == Nθ_points) ? 1 : j + 1
    jm = (j == 1) ? Nθ_points : j - 1

    radial =
        ((r[i+1]*ur[i+1,j,k]) -
         (r[i-1]*ur[i-1,j,k]))/(2Δr)

    angular =
        (uθ[i,jp,k] -
         uθ[i,jm,k])/(2Δθ)

    vertical =
        (uz[i,j,k+1] -
         uz[i,j,k-1])/(2Δz)

    return radial/ri +
           angular/ri +
           vertical

end

# ----------------------------------------------------------
# BUOYANCY FORCE
# ----------------------------------------------------------

function buoyancy(Tvalue)

    return g*β*(Tvalue - T₀)

end

println("Cylindrical finite-difference operators loaded successfully.")




 # ==========================================================
# BLOCK 3
# TENTATIVE VELOCITY SOLVER
# Computes ur*, uθ*, uz*
# ==========================================================
function tentative_velocity!()

println("Computing tentative velocities...")

# ----------------------------------------------------------
# Interior points
# ----------------------------------------------------------

for i = 2:Nr_points-1

    for j = 1:Nθ_points

        for k = 2:Nz_points-1

            # --------------------------------------------
            # Radial Momentum
            # --------------------------------------------

            adv_r =
                advection(ur,ur,uθ,uz,i,j,k)

            diff_r =
                ν*laplacian(ur,i,j,k)

            ri = max(r[i],Δr/2)

            centrifugal =
                uθ[i,j,k]^2/ri

            coupling =
                (2ν/ri^2)*
                dθ(uθ,i,j,k)

            ur_star[i,j,k] =
                ur[i,j,k]
                +
                dt*(
                -adv_r
                +
                diff_r
                +
                centrifugal
                -
                ν*ur[i,j,k]/ri^2
                -
                coupling
                )

            # --------------------------------------------
            # Azimuthal Momentum
            # --------------------------------------------

            adv_θ =
                advection(uθ,ur,uθ,uz,i,j,k)

            diff_θ =
                ν*laplacian(uθ,i,j,k)

            curvature =
                ur[i,j,k]*uθ[i,j,k]/ri

            coupling2 =
                (2ν/ri^2)*
                dθ(ur,i,j,k)

            uθ_star[i,j,k] =
                uθ[i,j,k]
                +
                dt*(
                -adv_θ
                +
                diff_θ
                -
                curvature
                -
                ν*uθ[i,j,k]/ri^2
                +
                coupling2
                )

            # --------------------------------------------
            # Axial Momentum
            # --------------------------------------------

            adv_z =
                advection(uz,ur,uθ,uz,i,j,k)

            diff_z =
                ν*laplacian(uz,i,j,k)

            buoy =
                buoyancy(T[i,j,k])

            uz_star[i,j,k] =
                uz[i,j,k]
                +
                dt*(
                -adv_z
                +
                diff_z
                +
                buoy
                )

        end

    end

end

println("Tentative velocities computed.")
return nothing
end




 # ==========================================================
# BLOCK 4
# CYLINDRICAL PRESSURE POISSON SOLVER
# ==========================================================
function pressure_poisson!()

println("Solving cylindrical Pressure Poisson equation...")

pressure_iterations = 20

for iter = 1:pressure_iterations

    p_old = copy(p)

    for i = 2:Nr_points-1

        ri = max(r[i], Δr/2)

        for j = 1:Nθ_points

            jp = (j == Nθ_points) ? 1 : j+1
            jm = (j == 1) ? Nθ_points : j-1

            for k = 2:Nz_points-1

                rhs =
                (ρ₀/dt) *
                divergence(
                    ur_star,
                    uθ_star,
                    uz_star,
                    i,j,k
                )

                radial =
                (
                    (p_old[i+1,j,k]-2p_old[i,j,k]+p_old[i-1,j,k])/(Δr^2)
                )
                +
                (
                    (p_old[i+1,j,k]-p_old[i-1,j,k])/
                    (2ri*Δr)
                )

                angular =
                (
                    p_old[i,jp,k]
                    -
                    2p_old[i,j,k]
                    +
                    p_old[i,jm,k]
                )/(ri^2*Δθ^2)

                vertical =
                (
                    p_old[i,j,k+1]
                    -
                    2p_old[i,j,k]
                    +
                    p_old[i,j,k-1]
                )/(Δz^2)

                residual =
                radial +
                angular +
                vertical -
                rhs

                diagonal =
                2/(Δr^2)
                +
                2/(ri^2*Δθ^2)
                +
                2/(Δz^2)

                p[i,j,k] =
                p_old[i,j,k]
                -
                residual/diagonal

            end
        end
    end

    # ------------------------------------------------------
    # Pressure Boundary Conditions
    # ------------------------------------------------------

    p[1,:,:] .= p[2,:,:]

    p[end,:,:] .= p[end-1,:,:]

    p[:,:,1] .= p[:,:,2]

    p[:,:,end] .= p[:,:,end-1]

end

println("Pressure solve complete.")

return nothing
end





# ==========================================================
# BLOCK 5
# VELOCITY PROJECTION
# Correct velocities using the pressure field
# ==========================================================
function project_velocity!()

println("Projecting velocity field...")

for i = 2:Nr_points-1

    ri = max(r[i], Δr/2)

    for j = 1:Nθ_points

        jp = (j == Nθ_points) ? 1 : j + 1
        jm = (j == 1) ? Nθ_points : j - 1

        for k = 2:Nz_points-1

            # --------------------------------------------------
            # Pressure gradients
            # --------------------------------------------------

            dpdr =
            (p[i+1,j,k] - p[i-1,j,k])/(2Δr)

            dpdθ =
            (p[i,jp,k] - p[i,jm,k])/(2Δθ)

            dpdz =
            (p[i,j,k+1] - p[i,j,k-1])/(2Δz)

            # --------------------------------------------------
            # Correct velocities
            # --------------------------------------------------

            ur[i,j,k] =
            ur_star[i,j,k]
            -
            dt/ρ₀ * dpdr

            uθ[i,j,k] =
            uθ_star[i,j,k]
            -
            dt/(ρ₀*ri) * dpdθ

            uz[i,j,k] =
            uz_star[i,j,k]
            -
            dt/ρ₀ * dpdz

        end
    end
end

println("Velocity projection complete.")

return nothing
end






# ==========================================================
# BLOCK 6
# TEMPERATURE UPDATE
# ==========================================================
function update_temperature!()

println("Updating temperature field...")

T_new = copy(T)

for i = 2:Nr_points-1

    for j = 1:Nθ_points

        for k = 2:Nz_points-1

            # ------------------------------------------
            # Advection
            # ------------------------------------------

            advT =
            advection(
                T,
                ur,
                uθ,
                uz,
                i,j,k
            )

            # ------------------------------------------
            # Diffusion
            # ------------------------------------------

            diffT =
            κ *
            laplacian(
                T,
                i,j,k
            )

            # ------------------------------------------
            # Explicit Euler update
            # ------------------------------------------

            T_new[i,j,k] =
            T[i,j,k]
            +
            dt*(
                -advT
                +
                diffT
            )

        end

    end

end

T .= T_new

println("Temperature updated.")

return nothing
end



# ==========================================================
# BLOCK 7
# CONCENTRATION UPDATE
# ==========================================================
function update_concentration!()

println("Updating concentration field...")

C_new = copy(C)

for i = 2:Nr_points-1

    for j = 1:Nθ_points

        for k = 2:Nz_points-1

            # ------------------------------------------
            # Advection
            # ------------------------------------------

            advC =
            advection(
                C,
                ur,
                uθ,
                uz,
                i,j,k
            )

            # ------------------------------------------
            # Diffusion
            # ------------------------------------------

            diffC =
            D *
            laplacian(
                C,
                i,j,k
            )

            # ------------------------------------------
            # Explicit Euler update
            # ------------------------------------------

            C_new[i,j,k] =
            C[i,j,k]
            +
            dt*(
                -advC
                +
                diffC
            )

        end

    end

end

C .= C_new

println("Concentration updated.")


return nothing
end





# ==========================================================
# BLOCK 8
# BOUNDARY CONDITIONS
# ==========================================================

function apply_boundary_conditions!()

    println("Applying boundary conditions...")

    # ----------------------------------------------------------
    # VELOCITY
    # ----------------------------------------------------------
    # Cylindrical wall: r = R
    # No-slip condition
    # ----------------------------------------------------------

    ur[end,:,:] .= 0.0
    uθ[end,:,:] .= 0.0
    uz[end,:,:] .= 0.0

    # ----------------------------------------------------------
    # Axis: r = 0
    # Regularity conditions
    # ----------------------------------------------------------

    ur[1,:,:] .= 0.0
    uθ[1,:,:] .= 0.0
    uz[1,:,:] .= uz[2,:,:]

    # ----------------------------------------------------------
    # Bottom plate: z = 0
    # No-slip condition
    # ----------------------------------------------------------

    ur[:,:,1] .= 0.0
    uθ[:,:,1] .= 0.0
    uz[:,:,1] .= 0.0

    # ----------------------------------------------------------
    # Top plate: z = H
    # No-slip condition
    # ----------------------------------------------------------

    ur[:,:,end] .= 0.0
    uθ[:,:,end] .= 0.0
    uz[:,:,end] .= 0.0

    # ----------------------------------------------------------
    # TEMPERATURE
    # ----------------------------------------------------------

    # Bottom plate: z = 0
    T[:,:,1] .= Tb

    # Top plate: z = H
    T[:,:,end] .= Tt

    # Cylindrical wall: r = R
    # Insulated wall: ∂T/∂r = 0
    T[end,:,:] .= T[end-1,:,:]

    # Axis: r = 0
    # Symmetry: ∂T/∂r = 0
    T[1,:,:] .= T[2,:,:]

    # ----------------------------------------------------------
    # CONCENTRATION
    # ----------------------------------------------------------
    # Zero-flux (Neumann) conditions
    # ----------------------------------------------------------

    # Axis: r = 0
    C[1,:,:] .= C[2,:,:]

    # Cylindrical wall: r = R
    C[end,:,:] .= C[end-1,:,:]

    # Bottom: z = 0
    C[:,:,1] .= C[:,:,2]

    # Top: z = H
    C[:,:,end] .= C[:,:,end-1]

    # ----------------------------------------------------------
    # PRESSURE
    # ----------------------------------------------------------
    # Homogeneous Neumann conditions
    # ----------------------------------------------------------

    p[1,:,:] .= p[2,:,:]
    p[end,:,:] .= p[end-1,:,:]

    p[:,:,1] .= p[:,:,2]
    p[:,:,end] .= p[:,:,end-1]

    println("Boundary conditions applied.")

    return nothing
end






# ==========================================================
# BLOCK 9
# TIME-STEPPING LOOP
# ==========================================================

println("Starting time integration...")

for n = 1:Nt

    tentative_velocity!()

    pressure_poisson!()

    project_velocity!()

    update_temperature!()

    update_concentration!()

    apply_boundary_conditions!()

    if n % 50 == 0
        println("Completed time step $n of $Nt")
    end

end

println("Time integration complete.")




# ==========================================================
# BLOCK 10
# TEMPERATURE HEAT MAP
# CYLINDRICAL COORDINATES
# ==========================================================

using Plots

println("Generating temperature heat map...")

# Select the middle azimuthal plane
j_mid = 1

# Temperature as a function of r and z
T_rz = T[:, j_mid, :]

heatmap(
    z,
    r,
    T_rz,
    xlabel = "z",
    ylabel = "r",
    title = "Temperature Distribution in the Cylinder",
    colorbar_title = "Temperature",
    aspect_ratio = :equal
)


# ==========================================================
# SAVE SIMULATION RESULTS
# ==========================================================

println("Saving simulation results...")

open("results/temperature.txt", "w") do io
    for k = 1:Nz_points
        for i = 1:Nr_points
            for j = 1:Nθ_points
                println(io, r[i], " ", θ[j], " ", z[k], " ", T[i,j,k])
            end
        end
    end
end

println("Temperature data saved to results/temperature.txt")
