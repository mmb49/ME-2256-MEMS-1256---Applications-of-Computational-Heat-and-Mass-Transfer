clear all
close all
clc

% Number dx
N = 10;

% Spatial domain definition
L = 10;
x = linspace(0,L,N);

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
    T_old(1) = 0;
    T_new(end) = 100;
    
    for i = 2:N-1
        T_new(i) = (T_old(i+1) + T_old(i-1))/2;
    end
    % Calculating the relative error between iterations
    error = max(T_new - T_old);
    % Resetting our temperature
    T_old = T_new;
    % Updating our loop counter
    counter = counter + 1;
end

% Plotting
plot(x,T_new,'-.b')
xlabel('x [length]')
ylabel('T [Temperature]')