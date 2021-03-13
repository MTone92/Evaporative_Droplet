module MyFunctions
using SpecialFunctions
using Roots
export LegendreP, GegenC, M1, M2, N1, N2, J, Jα, k1, k2

function LegendreP(τ,x)
    da = π/10000;
    s = 0;
    for a = 0.5*da:da:(π-0.5*da)
        s = s + da*(x+sqrt(Complex(x^2-1))*cos(a))^(-0.5-1im*τ);
    end
    
    return real(1/π*s);
end

function GegenC(τ,x)
    da = π/10000;
    s = 0;
    for a = (0.5*da):da:(π-0.5*da)
        s = s + da*((x+(x^2-1+0im)^0.5*cos(a))^2-1)/(x+(x^2-1+0im)^0.5*cos(a))^(1.5+1im*τ);
    end
    return real(1/(2*π*τ*1im)*s)
end

function M1(τ,βc)
    M1 = sin(βc)*sinh(τ*βc);
    return M1;
end

function M2(τ,βc)
    M2 = (τ^2-1)*sin(βc)*sinh(τ*βc)+2*τ*cos(βc)*cosh(τ*βc);
    return M2;
end

function N1(τ,βc)
    N1 = cos(βc)*sinh(τ*βc)-τ*sin(βc)*cosh(τ*βc);
    return N1;
end

function N2(τ,βc)
    N2 = (τ^2+1)*(cos(βc)*sinh(τ*βc)+τ*sin(βc)*cosh(τ*βc));
    return N2;
end

function J(α,βc)
    dτ = 0.01
    τ_max = 20
    coshα = cosh(α)
    J1 = 0
    for τ = dτ/2 : dτ : (τ_max - dτ/2)
        J1 = J1 + cosh(βc*τ)/cosh(π*τ)*tanh((π-βc)*τ)*LegendreP(τ,cosh(α))*τ
    end
    J = sin(βc)/2 + √2*(cosh(α)+cos(βc))^(3/2) * J1*dτ
    return J
end

function Jα(α,βc)
    dτ = 0.01
    τ_max = 20
    coshα = cosh(α)
    J1 = 0
    for τ = dτ/2 : dτ : (τ_max - dτ/2)
        J1 = J1 + cosh(βc*τ)/cosh(π*τ)*tanh((π-βc)*τ)*(τ^2+1/4)/sinh(α)*GegenC(τ,cosh(α))*τ
    end
    J = √2*(cosh(α)+cos(βc))^(3/2) * J1*dτ
    return J
end

function k1(τ,βc)
    k1 = (N2(τ,βc)*K(βc,τ)+N1(τ,βc)*K_bar(βc,τ))/(N2(τ,βc)*M1(τ,βc)+N1(τ,βc)*M2(τ,βc));
    return k1;
end

function k2(τ,βc)
    k2 = (M2(τ,βc)*K(βc,τ)-M1(τ,βc)*K_bar(βc,τ))/(N2(τ,βc)*M1(τ,βc)+N1(τ,βc)*M2(τ,βc));
    return k2;
end


end