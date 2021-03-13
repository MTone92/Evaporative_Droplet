function f = GegenC_myself(tau,x)
    da = pi/10000;
    s = 0;
    for a = (0.5*da):da:(pi-0.5*da)
        s = s + da*((x+(x.^2-1).^0.5*cos(a)).^2-1)./(x+(x.^2-1).^0.5*cos(a)).^(1.5+1i*tau);
    end
    f = real(1/(2*pi*tau*1i)*s);
end