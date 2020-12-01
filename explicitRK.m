
function [t, u] = explicitRK(func, y0, t0, tn, dt)

%Initial Conditions
n = round((tn - t0)/dt) + 1;
t = transpose(linspace(t0, tn, n));
u = zeros(n, 1);
u(1, 1) = y0;

%Define Butcher Array of desired RK Method: for example  the RK4 method
%has:
c = [0; 1/2; 1/2; 1]; %coefficients of Kj in second term of each Ki
b = [1/6 1/3 1/3 1/6]; %coefficients of Ki's*2 in un+1 term
a = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0]; %coefficients of dt in first term of Ki's

%Gives:
%[0  :                ]
%[1/2: 1/2            ]
%[1/2: 0   1/2        ]
%[1  : 0   0   1      ]
%[ - : 1/6 1/3 1/3 1/6]

%Note: one can change butcher array so long as the coeffients satisfy the
%conditions that were used to prove the method of order n that you are
%using.

%The following for loop effectively calculates un+1 = un + dt(sum(bi*ki))
for i = 1:1:n-1
    
    ks = zeros(length(b), 1);
    tn = t(i, 1);
    un = u(i, 1);
    
    %This for loop calculates: ki = f(tn+ci*dt, un + dt(sum(aij*kj))
    for k = 1:1:length(b)
        
        ci = c(k, 1);
        
        ks(k, 1) = func((tn + ci), (un + dt*dot(a(k,:), ks)));
    end
        
    u(i+1, 1) = u(i, 1) + dt*dot(b, ks);
end
