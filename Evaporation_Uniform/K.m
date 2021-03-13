function f = K(thetac,tau)
    K = -1./(sqrt(2).*cosh(pi.*tau)).*(2.*tau.*cot(thetac/2).*sinh(thetac.*tau)+cosh(thetac.*tau));
    f = K;
end