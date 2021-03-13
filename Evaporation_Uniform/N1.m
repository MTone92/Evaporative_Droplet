function f = N1(tau,thetac)
    N1 = cos(thetac).*sinh(tau.*thetac)-tau.*sin(thetac).*cosh(tau.*thetac);
    f = N1;
end