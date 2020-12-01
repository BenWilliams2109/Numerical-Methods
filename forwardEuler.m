%Forward Euler Method%
%clear
%clc
%close('all')

function [t, u] = forwardEuler(func, yInit, tInit, tFinal, tChange)
%figure how how many time instants t0, t1, ..., tn
n = round((tFinal-tInit)/tChange) + 1;
  
%instants t0, t1, ..., tn
t = linspace(tInit, tFinal, n);
    
%Initialize solution vector
u = zeros(n, 1);
    
%Initial value
u(1, 1) = yInit;

for i = 1:1:n-1
    u(i+1, 1) = u(i, 1) + (tChange * func(t(1, i), u(i, 1)));
end

%solvedinho = [u t];
%disp(solvedinho)
%plot(u1, t1, u2, t2, u3, t3)
%legend('dt=0.01', 'dt=0.25', 'dt=0.5')
%title('Forward Euler Approximation')
%ylabel('Approximate Solution')
%xlabel('Time')
end

%figure(1);
%plot(t1,u1,t2,u2,t3,u3);
%title(’Forward Euler method’);
%legend(’dt=0.01’,’dt=0.25’,’dt=0.5’)
%figure(2);
%plot(t4,u4,t5,u5,t6,u6);
%title(’Backward Euler method’);
%legend(’dt=0.01’,’dt=0.25’,’dt=0.5’)
