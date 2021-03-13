function f = k2(thetac,tau)
    k2 = (M2(tau,thetac).*K(thetac,tau)-M1(tau,thetac).*K_bar(thetac,tau))./(N2(tau,thetac).*M1(tau,thetac)+N1(tau,thetac).*M2(tau,thetac));
    f = k2;
end