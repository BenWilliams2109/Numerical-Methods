
function [t, u] = CrankNicholson(func, y0, t0, tn, dt)
    
n = round((tn - t0)/dt) + 1;
t = transpose(linspace(t0, tn, n));
u = zeros(n, 1);
u(1, 1) = y0;

options = optimset;
options.Display = 'off';
options.TolFun = 1.e-12;
options.MaxFunEvals = 10000;

for i = 1:1:n-1
    
    tnplus1 = t(i+1, 1);
    un = u(i, 1);
    
    %fzero takes arguement unplus1 and solves rearranged backward euler
    %method, such that it is all rearranged to equal 0, for unplus1.
    %Starting at un and taking options.
    %Replace equation with "backwardEulerFunc(arguaments as shown below)" to
    %spit out vector z containing error at each step.
    unplus1 = fzero(@(unplus1) (unplus1 - un - (dt/2)*(func(tn, un) + func(tnplus1, unplus1))), un, options);
    
    u(i+1, 1) = unplus1;
    
end
end