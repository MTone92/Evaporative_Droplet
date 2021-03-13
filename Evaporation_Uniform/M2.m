function f = M2(tau,thetac)
    M2 = (tau.^2-1).*sin(thetac).*sinh(tau.*thetac)+2*tau.*cos(thetac).*cosh(tau.*thetac);
    f = M2;
end