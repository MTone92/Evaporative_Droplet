function f = k1(thetac,tau)
    k1 = (N2(tau,thetac).*K(thetac,tau)+N1(tau,thetac).*K_bar(thetac,tau))./(N2(tau,thetac).*M1(tau,thetac)+N1(tau,thetac).*M2(tau,thetac));
    f = k1;
end