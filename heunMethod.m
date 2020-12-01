function [t, u] = heunMethod(func, y0, t0, tn, dt)

n = round((tn - t0)/dt) + 1;
t = transpose(linspace(t0, tn, n));
u = zeros(n, 1);
u(1, 1) = y0;

for i = 1:1:n-1
    
    tn = t(i, 1);
    tnplus1 = t(i+1, 1);
    un = u(i, 1);

    unplus1 = un + (dt/2)*(func(tn, un) + func(tnplus1, (un + dt*func(tn, un))));
    u(i+1, 1) = unplus1;
    
end
end