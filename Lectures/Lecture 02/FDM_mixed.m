clear all
close all
clc

% Number dx
N = 5;

% Spatial domain definition
L = 1;
x = linspace(0,L,N);
dx = L/(N-1);
lambda = 1;
A = 1;

% Temperature array preallocation
T_old = zeros(1,N);
T_new = zeros(1,N);

% Defining relative error
error = 1;
% Defining a loop counter
counter = 0;

% Solving for interior temperatures
while error >= 1e-5
    % Bounday conditions
    q_dp = 100;
    T_new(end) = 300;
    T_0 = T_old(2) + (2*dx*q_dp)/(lambda*A);
    for i = 1:N-1
        if i == 1
           T_new(i) = (T_old(i+1) + T_0)/2; 
        else
           T_new(i) = (T_old(i+1) + T_old(i-1))/2;
        end
    end
    % Calculating the relative error between iterations
    error = max(T_new - T_old);
    % Resetting our temperature
    T_old = T_new;
    % Updating our loop counter
    counter = counter + 1;
end

x_analytic = linspace(0,L,100);
T_analytic = @(x) 300 + (q_dp/lambda)*(L - x); 

% Plotting
plot(x,T_new,'*b')
hold on
plot(x_analytic,T_analytic(x_analytic),'k')
xlabel('x [length]')
ylabel('T [Temperature]')