% This case is based on the following assumptions:
% 1. Spherical cap shape.
% 2. Pinned contact line.
% 3. Uniform evaperative flux.

% Tips:

% 1.Based on the experiment with tau = 100,200 and 300, the stream value is
% 1.1907e-24, 1.1925e-24, 1.1925e-24, respectively. Thus, we can just use
% tau = 200 as the upperbound of the interval.

% 2. The maximum interval of the integration of the SPI is tau = (0,226),
% it because that cosh(pi*tau) will exceed the maximum number can be stored
% in memory and become cosh(pi*tau) = Inf, which is not a number and caut
% error.

% The main program:
rho = 1;
R = 1;
J0 = 1;
%thetac = pi/2;

r_1 = linspace(0,R,101);
z_1 = sqrt(R^2/sin(thetac)^2-r_1.^2)-R*cot(thetac);
[r,z] = meshgrid(r_1);
for i = 1:1:101
    zmax = sqrt(R^2/sin(thetac)^2-r(1,i)^2)-R*cot(thetac);
    z(:,i) = linspace(0,zmax,101);
end
alpha_p = acoth((r.^2+z.^2+R.^2)./(2*R.*r));
theta_p = acos(R.*sinh(alpha_p)./r-cosh(alpha_p));
streamfun = zeros(101);
alpha = alpha_p(2:101,2:100);
theta = theta_p(2:101,2:100);

streamfun(2:101,2:100) = PSI(thetac,theta,alpha);
