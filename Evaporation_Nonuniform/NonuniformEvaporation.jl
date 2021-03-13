# This program is used to calculate the dimensionless Stream function inside a droplet by Nonuniform evaporation.

# Notes for coding:
# 1. Two important function in this code is P_(-1/2 + iτ) and C_(1/2 + iτ) ^(-1/2), which are LegendreP and GegenC.
# 2.1 Feature 1 of LegendreP(τ,cosh(α)): converge to 0 and vibrate less with α increase. When α = 0 (cosh(α)=1), LegendreP(τ,cosh(α)) reaches maximum value, 1.
# 2.2 Feature 2 of LegendreP(τ,cosh(α)): vibrate in range (-1, 1) constantly with τ increase. When α is larger, vibrate more frequently in τ-direction with lower value.
# 3.1 Feature 1 of GegenC(τ,cosh(α)): diverge and vibrate more frequently with α increase
# 3.2 Feture 2 of GegenC(τ,cosh(α)): converge to 0 and vibrate in constant frequency with τ increase
# 4. When calculate dβcdt, it is fine to integrate in α in range [0, 20], since J(α=20) will be very smalle due to exp(7-α) term.
# 5.1 When consider the integration in dβcdt, ψc and ψc_bar, sinh(α)/(cosh(α)+cos(βc))^2*J(α,βc) is a smooth curve, GegenC(τ,cosh(α)) vibrate with α increase.
#       However, the maximum frequency of GegenC(τ,cosh(α)) via α (@α=20) is 2. Thus, dα=0.01 is sufficient to offer accurate integral.
# 6. When finding J, flux at the droplet surface, it is suggested to integrate from 0 to 7 in α-direction.

# Main program:
# Import the necessary packages before running the main program:
using Base.Threads
using Distributed
addprocs()
@everywhere include("Module MyFunctions.jl")        # Include the file directory and file name for using the local modules
@everywhere using .MyFunctions                      # The dot before the module name is used to indicate it is a local module
@everywhere using SpecialFunctions
using Dates
using MAT

StartTime = Dates.now()
# Check the threads' number:
N_t = nthreads()
N_w = nworkers()
println("Number of Threads = ",N_t)
println("Number of workers = ",N_w)

# About the droplet:
βc = π/20                            # Contact angle
R = 1e-3                            # Droplet Radius
T = 293.15                          # Room Temperature
w = 14.62*10^(-3);                  # This is absolute humidity under Temp = 20C
ρ_w = 998.2                         # Desnity of water @ 20C
ρ_v = 1.194                         # Density of vaper-air mixture @ 20C
Humidity = 0                        # Relative humidity at far field

P_vapor = 2.339                     # Saturated vapor pressure under Temp = 20C.
P_total = 101.353                   # Atmospheric pressure at sea level.
w = 0.662*P_vapor/(P_total-P_vapor)
D = 22.5*10^(-6)*(T/273.15)^1.8     # Binary diffusion coefficient of vapor in air.
J0 = ρ_v*D*w*(1-Humidity)/R         # Reference flux at droplet surface
ψ0 = J0*R^2/ρ_w                     # Reference stream function value

# Define calculation grid:
r1 = [0:0.01:1;]
z1 = sqrt.(1/sin(βc)^2 .- r1.^2) .- cot(βc)
lengthr = length(r1)
lengthz = length(z1)
r = repeat(r1',outer=(lengthz,1))
z = repeat(z1,outer=(1,lengthr))
for i = 1:1:100
    zmax = sqrt(1/sin(βc)^2-r[1,i]^2) - cot(βc)
    z[:,i] = [0:zmax/100:zmax+zmax/200;]
end
z[:,101] = zeros(lengthz)
alpha_p = acoth.((r.^2 .+ z.^2 .+ 1)./(2r))
beta_p = sinh.(alpha_p)./r .- cosh.(alpha_p)
beta_p[1,:] = ones(1,101)
beta_p = acos.(beta_p)
α = alpha_p[2:101,2:100]
β = beta_p[2:101,2:100]
Ψ = zeros(lengthz,lengthr)

# Calculations:
# Caution: Need to be checked
# Find J(α,βc)
# dα1 = 1                         # For test
dα1 = 0.001
α1 = [dα1/2 : dα1 : 7-dα1/4;]
α2 = [7+dα1/2 : dα1 : 7-dα1/4;]
lengthα1 = length(α1)
J1 = zeros(size(α1))
J2 = J(7,βc).*exp.(7 .- α2)

Threads.@threads for i = 1:1:lengthα1
    J1[i] = J(α1[i],βc)
end
α3 = [α1; α2]
J3 = [J1; J2]
lengthα3 = length(α3)

# Find Vβ(α,βc)
dβcdt = @distributed (+) for i = 1:1:lengthα3
    sinh(α3[i])/(cosh(α3[i])+cos(βc))^2 * J3[i]
end
dβcdt = -2*(1+cos(βc))^2 * dβcdt * dα1

Vβc = dβcdt./(cosh.(α3).+cos(βc)) .+ J3

# Find ψβ(α,βc)
ψβc = zeros(lengthα3)
ψβc[1] = sinh(α3[1])/(cosh(α3[1])+cos(βc))^2 * Vβc[1]
for i = 2:1:lengthα3
    ψβc[i] = ψβc[i-1] + sinh(α3[i])/(cosh(α3[i])+cos(βc))^2 * Vβc[i]
end
ψβc = -ψβc * dα1

# Caution: Need to be checked
# Find dJdα(α,βc)
Jα1 = zeros(lengthα1)
Threads.@threads for i = 1:1:lengthα1
    Jα1[i] = Jα(α1[i],βc)
end

dJdα1 = zeros(lengthα1)
Threads.@threads for i = 1:1:lengthα1
    dJdα1[i] = (J1[i]-sin(βc)/2) * 3/2*sinh(α1[i])/(cosh(α1[i])+cos(βc)) + Jα1[i]
end

dJdα2 = -J2
dJdα3 = [dJdα1; dJdα2]

# Find ψc_bar(α,βc)
ψβc_bar = zeros(lengthα3)
Threads.@threads for i = 1:1:lengthα3
    ψβc_bar[i] = -(cosh(α3[i])+cos(βc))^(-1/2) * 
                    (sinh(α3[i])/(cosh(α3[i])+cos(βc)) * (sinh(α3[i])*J3[i] + (cosh(α3[i])+cos(βc))*dJdα3[i]) +
                    3/2*ψβc[i] * (cos(βc)*(cosh(α3[i])+cos(βc)) - sin(βc)^2/2))
end


# Calculate Psi:
# dτ1 = 0.1               # For test
# dτ2 = 1                 # For test
# dτ3 = 10                # For test
dτ1 = 0.001
dτ2 = 0.01
dτ3 = 0.1
Ψ1 = @distributed (+) for τ = dτ1/2 : dτ1 : 2-dτ1/4
# This loop is used to interate in τ-direction
    # Find K(τ,βc)
    K = 0
    for i = 1:1:lengthα3
        K += (cosh(α3[i])+cos(βc))^(3/2)/sinh(α3[i]) * ψβc[i] * GegenC(τ,cosh(α3[i]))
    end
    K = τ*(τ^2+1/4)*tanh(π*τ)*K*dα1

    # Find K_bar
    K_bar = 0
    for i = 1:1:lengthα3
        K_bar += ψβc_bar[i]/sinh(α3[i]) * GegenC(τ,cosh(α3[i]))
    end
    K_bar = τ*(τ^2+1/4)*tanh(π*τ)*K_bar*dα1

    # Find ψ1
    k1 = (N2(τ,βc)*K + N1(τ,βc)*K_bar)/(N2(τ,βc)*M1(τ,βc) + N1(τ,βc)*M2(τ,βc))
    k2 = (M2(τ,βc)*K - M1(τ,βc)*K_bar)/(N2(τ,βc)*M1(τ,βc) + N1(τ,βc)*M2(τ,βc))
    K_main = k1.*sin.(β).*sinh.(τ.*β) .+ k2.*(cos.(β).*sinh.(τ.*β) .- τ.*sin.(β).*cosh.(τ.*β))
    
    K_main.*GegenC.(τ,cosh.(α))
end
Ψ1 = Ψ1 .* dτ1

Ψ2 = @distributed (+) for τ = 2+dτ2/2 : dτ2 : 10-dτ2/4
# This loop is used to interate in τ-direction
    # Find K(τ,βc)
    K = 0
    for i = 1:1:lengthα3
        K += (cosh(α3[i])+cos(βc))^(3/2)/sinh(α3[i]) * ψβc[i] * GegenC(τ,cosh(α3[i]))
    end
    K = τ*(τ^2+1/4)*tanh(π*τ)*K*dα1

    # Find K_bar
    K_bar = 0
    for i = 1:1:lengthα3
        K_bar += ψβc_bar[i]/sinh(α3[i]) * GegenC(τ,cosh(α3[i]))
    end
    K_bar = τ*(τ^2+1/4)*tanh(π*τ)*K_bar*dα1

    # Find ψ1
    k1 = (N2(τ,βc)*K + N1(τ,βc)*K_bar)/(N2(τ,βc)*M1(τ,βc) + N1(τ,βc)*M2(τ,βc))
    k2 = (M2(τ,βc)*K - M1(τ,βc)*K_bar)/(N2(τ,βc)*M1(τ,βc) + N1(τ,βc)*M2(τ,βc))
    K_main = k1.*sin.(β).*sinh.(τ.*β) .+ k2.*(cos.(β).*sinh.(τ.*β) .- τ.*sin.(β).*cosh.(τ.*β))
    
    K_main.*GegenC.(τ,cosh.(α))
end
Ψ2 = Ψ2 .* dτ2

Ψ3 = @distributed (+) for τ = 10+dτ3/2 : dτ3 : 220-dτ3/4
# This loop is used to interate in τ-direction
    # Find K(τ,βc)
    K = 0
    for i = 1:1:lengthα3
        K += (cosh(α3[i])+cos(βc))^(3/2)/sinh(α3[i]) * ψβc[i] * GegenC(τ,cosh(α3[i]))
    end
    K = τ*(τ^2+1/4)*tanh(π*τ)*K*dα1

    # Find K_bar
    K_bar = 0
    for i = 1:1:lengthα3
        K_bar += ψβc_bar[i]/sinh(α3[i]) * GegenC(τ,cosh(α3[i]))
    end
    K_bar = τ*(τ^2+1/4)*tanh(π*τ)*K_bar*dα1

    # Find ψ1
    k1 = (N2(τ,βc)*K + N1(τ,βc)*K_bar)/(N2(τ,βc)*M1(τ,βc) + N1(τ,βc)*M2(τ,βc))
    k2 = (M2(τ,βc)*K - M1(τ,βc)*K_bar)/(N2(τ,βc)*M1(τ,βc) + N1(τ,βc)*M2(τ,βc))
    K_main = k1.*sin.(β).*sinh.(τ.*β) .+ k2.*(cos.(β).*sinh.(τ.*β) .- τ.*sin.(β).*cosh.(τ.*β))
    
    K_main.*GegenC.(τ,cosh.(α))
end
Ψ3 = Ψ3 .* dτ3

Ψ[2:101,2:100] = (cosh.(α) .+ cos.(β)).^(-3/2) .* (Ψ1 .+ Ψ2 .+ Ψ3)


# Find out the time elapse for running this whole program:
EndTime = Dates.now()
ΔTime = EndTime - StartTime
ΔTime = convert(Dates.Nanosecond,Dates.Millisecond(ΔTime))
ΔTime = Dates.Time(ΔTime)

# Write out other parameters:
file = open("Parameters.txt","w")
write("Parameters.txt",
        "Number of Threads = ","$N_t\n",
        "Number of workers = ","$N_w\n",
        "ΔTime = $ΔTime")
close(file)

matwrite("NonuniformEvaporation_b60_t200.mat",
    Dict("Streamfun" => Ψ, "r" => r, "z" => z, "alpha" => α, "beta" => β,
        "a3" => α3, "J3" => J3, "dJda3" => dJdα3, "V_bc" => Vβc,
        "psi_bc" => ψβc, "psi_bc_bar" => ψβc_bar, "beta" => βc,
        "Psi1" => Ψ1, "Psi2" => Ψ2, "Psi3" => Ψ3))