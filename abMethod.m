function [t,u] = abMethod(func, y0, t0, tn, dt)

n = round((tn - t0)/dt) + 1;
t = transpose(linspace(t0, tn, n));
u = zeros(n, 1);
u(1, 1) = y0;

%take p = 2, approximate f by the polynomial of degree 2 at tn, tn-1 and
%tn-2 to obtain the three-step method (AB3)

%Initialize u0, u1, u2 using RK method of similar order
[t(1:3), u(1:3)] = explicitRK(func, y0, t0, t(3,1), dt);


for i = 3:1:n-1
    
    %Formula
    u(i+1, 1) = u(i, 1) + (dt/12)*(23*func(t(i,1), u(i,1)) - ...
    16*func(t(i-1,1), u(i-1,1)) + 5*func(t(i-2,1), u(i-2,1)));
    
end
