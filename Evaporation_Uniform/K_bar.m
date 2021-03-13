function f = K_bar(thetac,tau)
    K_bar = 1./(sqrt(8).*cosh(pi.*tau).*(sin(thetac/2))^2).*(2*tau.*sinh(thetac.*tau).*cot(thetac/2).*(tau.^2.*(cos(thetac)-1)-2.*cos(thetac)-1)+cosh(thetac.*tau).*(tau.^2.*(5*cos(thetac)+7)-cos(thetac)+1));
    f = K_bar;
end