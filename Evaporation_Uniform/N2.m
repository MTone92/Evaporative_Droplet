function f = N2(tau,thetac)
    N2 = (tau.^2+1).*(cos(thetac)*sinh(tau.*thetac)+tau.*sin(thetac).*cosh(tau.*thetac));
    f = N2;
end