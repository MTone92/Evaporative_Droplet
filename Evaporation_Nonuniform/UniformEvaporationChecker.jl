# 90 degree checker:

βc = π/2
τ = 2.5
dα4 = 0.01                         # For test
α4 = [dα4/2:dα4:3-dα4/4;]
lengthα4 = length(α4)

# Value by uniform evaporation:
J4 = 1
dβcdt4 = -2*(1 + cos(βc))
Vβc4 = 1 .- 2(1 .+ cos(βc))./(cosh.(α4))
ψβc4 = (1 .- (1 + cos(βc))./(cosh.(α4) .+ cos(βc))) ./ (cosh.(α4) .+ cos(βc))
ψβc_bar4 = - .√(cosh.(α4) .+ cosh(βc)) .+ cos(βc)/2 ./ sqrt.(cosh.(α4) .+ cos(βc)) .+ 
            (4-(cos(βc)-3)^2/4)./(cosh.(α4) .+ cos(βc)).^(3/2) .-
            3sin(βc)^2*(1+cos(βc))/4 ./ (cosh.(α4) .+ cos(βc)).^(5/2)
K4 = -1 ./ (sqrt(2).*cosh(π.*τ)).*(2τ.*cot(βc/2).*sinh(βc.*τ) .+ cosh(βc.*τ))
K_bar4 = 1. /(sqrt(8).*cosh(π.*τ).*(sin(βc/2))^2).*(2τ.*sinh(βc.*τ).*cot(βc/2).*(τ.^2(cos(βc)-1)-2cos(βc)-1) +
            cosh(βc.*τ).*(τ.^2 .* (5*cos(βc)+7)-cos(βc)+1))

# Value by diffusive evaporation:
J5 = zeros(lengthα4)
Threads.@threads for i = 1:1:lengthα4
    J5[i] = J(α4[i],βc)
end

# Find Vβ(α,βc)
dβcdt5 = @distributed (+) for i = 1:1:lengthα4
    sinh(α4[i])/(cosh(α4[i])+cos(βc))^2 * J5[i]
end
dβcdt5 = -2*(1+cos(βc))^2 * dβcdt5 * dα4

Vβc5 = dβcdt5./(cosh.(α4).+cos(βc)) .+ J5

# Find ψβ(α,βc)
ψβc5 = zeros(lengthα4)
ψβc5[1] = sinh(α4[1])/(cosh(α4[1])+cos(βc))^2 * Vβc5[1]
for i = 2:1:lengthα4
    ψβc5[i] = ψβc5[i-1] + sinh(α4[i])/(cosh(α4[i])+cos(βc))^2 * Vβc5[i]
end
ψβc5 = -ψβc5 * dα4

# Caution: Need to be checked
# Find dJdα(α,βc)
Jα5 = zeros(lengthα4)
Threads.@threads for i = 1:1:lengthα4
    Jα5[i] = Jα(α4[i],βc)
end

dJdα5 = zeros(lengthα4)
Threads.@threads for i = 1:1:lengthα4
    dJdα5[i] = (J5[i]-sin(βc)/2) * 3/2*sinh(α4[i])/(cosh(α4[i])+cos(βc)) + Jα5[i]
end

# Find ψc_bar(α,βc)
ψβc_bar5 = zeros(lengthα4)
Threads.@threads for i = 1:1:lengthα4
    ψβc_bar5[i] = -(cosh(α4[i])+cos(βc))^(-1/2) * 
                    (sinh(α4[i])/(cosh(α4[i])+cos(βc)) * (sinh(α4[i])*J5[i] + (cosh(α4[i])+cos(βc))*dJdα5[i]) +
                    3/2*ψβc5[i] * (cos(βc)*(cosh(α4[i])+cos(βc)) - sin(βc)^2/2))
end

K5 = 0
for i = 1:1:lengthα4
    K5 += (cosh(α4[i])+cos(βc))^(3/2)/sinh(α4[i]) * ψβc5[i] * GegenC(τ,cosh(α4[i]))
end
K5 = τ*(τ^2+1/4)*tanh(π*τ)*K5*dα4

K_bar5 = 0
for i = 1:1:lengthα4
    K_bar5 += ψβc_bar4[i]/sinh(α4[i]) * GegenC(τ,cosh(α4[i]))
end
K_bar5 = τ*(τ^2+1/4)*tanh(π*τ)*K_bar5*dα4

# Test from real calculation:
K = @distributed (+) for i = 1:1:lengthα3
    (cosh(α3[i])+cos(βc))^(3/2)/sinh(α3[i]) * ψβc[i] * GegenC(τ,cosh(α3[i]))
end
K = τ*(τ^2+1/4)*tanh(π*τ)*K*dα1

K_bar = @distributed (+) for i = 1:1:lengthα3
    ψβc_bar[i]/sinh(α3[i]) * GegenC(τ,cosh(α3[i]))
end
K_bar = τ*(τ^2+1/4)*tanh(π*τ)*K_bar*dα1

KK1 = zeros(lengthα3)
Threads.@threads for i = 1:lengthα3
    KK1[i] = (cosh(α3[i])+cos(βc))^(3/2)/sinh(α3[i]) * ψβc[i] * GegenC(τ,cosh(α3[i]))
end

KK_bar1 = zeros(lengthα3)
Threads.@threads for i = 1:lengthα3
    KK_bar1[i] = ψβc_bar[i]/sinh(α3[i]) * GegenC(τ,cosh(α3[i]))
end

JJ3 = [1;J3];
dJJdα3 = J3;
dJJdα3[1] = (JJ3[2] - JJ3[1])/dα1;
dJJdα3[end] = (JJ3[end] - JJ3[end-1])/dα1;
Threads.@threads for i = 2:lengthα3-1
    dJJdα3[i] = (JJ3[i+1] - JJ3[i-1])/ 2dα1;
end