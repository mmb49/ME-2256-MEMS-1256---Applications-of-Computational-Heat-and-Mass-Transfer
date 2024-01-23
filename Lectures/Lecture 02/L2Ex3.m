% Lecture 2 Example #3
clear all
close all
clc

% Defining our constants:
L = 0.01; % [m]
Gamma = @(T) (5.238086549608868e-17).*T.^6 + ...
    (-2.927636770231909e-13).*T.^5 + ...
    (5.844390241944433e-10).*T.^4 + ...
    (-5.642804450717544e-7).*T.^3 + ...
    (0.0002909446395983974).*T.^2 + ...
    (-0.08063418038142083).*T + (11.00293123390308); % [W/m-K]
A_e = 2.5e-5; % [m^2]
A_w = A_e;

% Defining our domain:
N = 5; % number of control volumes
n = N+2; % number of nodes
delta_x_w = L/N;
delta_x_e = L/N;
dx = L/N;

% Defining our variables of interest:
T = 300*ones(1,(N+2));
coeff = zeros(N,N);
b = zeros(N,1);

% Defining our boundary conditions:
T(1) = 300; % [K]
T(end) = 650; % [K]

a_W = zeros(1,N);
a_E = zeros(1,N);
a_P = zeros(1,N);

% Defining our error criteria
error = 1;
counter = 0;

while error >= 1e-5
    % Defining our coefficients for the interior of our domain (CVs 1:N):
    for i = 1:N-1
        a_E(i) = A_e*(Gamma(T(i)) + Gamma(T(i+1)))/(2*delta_x_e);
    end
    a_E(end) = A_e*(2*Gamma(T(end)))/(delta_x_e);
    a_W(1) = A_w*(2*Gamma(T(1)))/(delta_x_w);
    for i = 2:N
        a_W(i) = A_w*(Gamma(T(i-1)) + Gamma(T(i)))/(2*delta_x_w);
    end
    for i = 1:N
        a_P(i) = a_E(i) + a_W(i);
    end
    
    % Populating our coefficient matrix
    for i = 1:N
        for j = 1:N
            if i == j
                coeff(i,j) = a_P(i);
            end
            if i + 1 == j
                coeff(i,j) = -a_E(i);
            end
            if j + 1 == i
                coeff(i,j) = -a_W(i);
            end
        end
    end
    
    % Next we populate the RHS matrix
    b(1) = a_W(1)*T(1);
    b(end) = a_E(end)*T(end);
    
    % Obtaining the solution through matrix inversion
    temp = coeff\b;
    
    % Defining maximum relative error
    error = abs(max(temp - transpose(T(2:end-1))));
    
    % Redefining our temperatures
    T(2:end-1) = temp;
    
    counter = counter + 1;
end

fprintf('T(0) = %.1f [K]\n',T(1))
for i = 1:N
    fprintf('T(%i) = %.1f [K]\n',i,temp(i))
end
fprintf('T(%i) = %.1f [K]\n',(N+1),T(end))

% Defining nodal points - where temperature is solved (used for plotting)
nodes = zeros(1,n);
nodes(1) = 0;
nodes(2) = dx/2;
for i = 3:N
    nodes(i) = nodes(i-1) + dx;
end
nodes(end-1) = L - (dx/2);
nodes(end) = L;

plot(nodes,[T(1);temp(:);T(end)],'-*b')
ylabel('Temperature [K]')
xlabel('Distance [m]')

