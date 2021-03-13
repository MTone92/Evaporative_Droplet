function f = GegenC(tau,x)
    fun = @(a) ((x+(x^2-1)^0.5*cos(a)).^2-1)./(x+(x^2-1)^0.5*cos(a)).^(1.5+1i.*tau);
    s = integral(fun,0,pi);
    f = real(1/(2*pi.*tau*1i).*s);
end