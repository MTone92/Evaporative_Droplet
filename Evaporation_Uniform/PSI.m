function f = PSI(thetac,theta,alpha)
    s = 0;
    dtau = 0.1;
    nodes = 30;
    spmd(nodes)
    for tau = ((labindex-1)*7+dtau/2):dtau:(labindex*7-dtau/2)
        s = s + dtau*(K_main(thetac,tau,theta).*GegenC_myself(tau,cosh(alpha)));
    end
    end

    g = 0;
    for i = 1:1:nodes
        g = g + s{i};
    end

    f = (cosh(alpha)+cos(beta)).^(-1.5).*g;
end