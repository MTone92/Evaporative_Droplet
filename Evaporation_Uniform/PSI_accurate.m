function f = PSI_accurate(thetac,theta,alpha)
    s = 0;
    dtau1 = 0.001;
    dtau2 = 0.01;
    dtau = 0.1;
    for tau = 0.0005:0.001:1.9995
        s = s + dtau1*(K_main(thetac,tau,theta).*(cosh(alpha)+cos(theta))^(-1.5).*GegenC(tau,cosh(alpha)));
    end
    
    for tau = 2.005:0.01:9.995
        s = s + dtau2*(K_main(thetac,tau,theta).*(cosh(alpha)+cos(theta))^(-1.5).*GegenC(tau,cosh(alpha)));
    end
    
    for tau = 10.05:0.1:220.05
        s = s + dtau*(K_main(thetac,tau,theta).*(cosh(alpha)+cos(theta))^(-1.5).*GegenC(tau,cosh(alpha)));
    end
    f = s;
end