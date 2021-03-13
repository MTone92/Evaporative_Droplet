function f = K_main(thetac,tau,theta)
    K_main = k1(thetac,tau).*sin(theta).*sinh(tau.*theta)+k2(thetac,tau).*(cos(theta).*sinh(tau.*theta)-tau.*sin(theta).*cosh(tau.*theta));
    f = K_main;
end