% Lecture 3 Example #1
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant dx

% Defining our constants:
L = 1; % [m]
Gamma_e = 1; % [W/m-K]
Gamma_w = 1; % [W/m-K]
A_e = 1; % [m^2]
A_w = A_e;

% Defining our domain:
N = 100; % number of control volumes
n = N+2; % number of nodes
delta_x_w = L/N;
delta_x_e = L/N;
dx = L/N;

% Defining nodal points - where temperature is solved (used for plotting)
nodes = zeros(1,n);
nodes(1) = 0;
nodes(2) = dx/2;
for i = 3:N
    nodes(i) = nodes(i-1) + dx;
end
nodes(end-1) = L - (dx/2);
nodes(end) = L;
    
% Defining our variables of interest:
T = zeros(1,n);
coeff = zeros(N,N);
b = zeros(N,1);

% Defining our boundary conditions:
T(1) = 100; % [K]
T(end) = 200; % [K]

% Defining our coefficients for the interior of our domain (CVs 2:N-1):
a_W = A_w*(Gamma_w/delta_x_w);
a_E = A_e*(Gamma_e/delta_x_e);
a_P = a_E + a_W;

% Populating our coefficient matrix. We first start with the interior
% diagonal values (i.e. a_P for CVs 2:N-1

for i = 2:N-1
    for j = 2:N-1
        if i == j
            coeff(i,j) = a_P; 
        end
    end
end

% Next we have to populate a_P for CVs 1 and N. Note for CV 1 we use the
% modified a_W expression, which considers delta_x_w/2
coeff(1,1) = a_E + (A_w*((2*Gamma_w)/dx));
% Note for CV 2 we use the modified expression for a_E, which considers
% delta_x_e/2
coeff(end,end) = a_W + (A_e*((2*Gamma_e)/dx));

% Next we populate the off-diagonals, i.e. a_E and a_W
for i = 1:N
    for j = 1:N
        if i + 1 == j
           coeff(i,j) = -a_E; 
        end
        if j + 1 == i
            coeff(i,j) = -a_W;
        end
    end
end

% Next we populate the RHS matrix
b(1) = (A_w*((2*Gamma_w)/dx))*T(1);
b(end) = (A_e*((2*Gamma_e)/dx))*T(end);

% Obtaining the solution through matrix inversion
temp = coeff\b;

temperatures = horzcat(T(1),transpose(temp(:)),T(end));

fprintf('Using a uniform grid:\n')
fprintf('T(0) = %.1f [K] at x=%.3f\n',T(1),nodes(1))
for i = 1:N
   fprintf('T(%i) = %.1f [K] at x=%.3f\n',i,temp(i),nodes(i+1))
end
fprintf('T(%i) = %.1f [K] at x=%.3f\n\n',(N+1),T(end),nodes(end))

plot(nodes,temperatures,'-ok')
hold on
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable dx

clearvars delta_x_w delta_x_e a_W a_E a_P coeff b T temp temperatures dx

% Defining our domain:
n = N+2; % number of nodes

% Defining nodal points
x = linspace(0,L,N);
for i = 1:N
    dx(i) = (L/N)*x(i) + (L/(2*N));
end
% Continuity check
if sum(dx(:)) - 1 >= 10^-15
    fprintf('domain discretization error.')
end

% Defining nodal points - where temperature is solved (used for plotting
% and defining delta_x_w and delta_x_e
nodes(1) = 0;
nodes(2) = dx(1)/2;
for i = 3:N
    nodes(i) = nodes(i-1) + (dx(i-2) + dx(i-1))/2;
end
nodes(end-1) = L - (dx(end)/2);
nodes(end) = L;

% Defining the nodal spacing, delta_x_w and delta_x_e
for i = 1:N
    delta_x_w(i) = nodes(i+1) - nodes(i);
    delta_x_e(i) = nodes(i+2) - nodes(i+1);
end

% Defining our variables of interest:
T = zeros(1,n);
coeff = zeros(N,N);
b = zeros(N,1);

% Defining our boundary conditions:
T(1) = 100; % [K]
T(end) = 200; % [K]

% Defining our coefficients a_W, a_P, and a_E
for i = 1:N
   a_E(i) = A_e*(Gamma_e/delta_x_e(i));
   a_W(i) = A_w*(Gamma_w/delta_x_w(i));
   a_P(i) = a_E(i) + a_W(i);
end

% Populating our coefficient matrix. We first start with the interior
% diagonal values (i.e. a_P for CVs 2:N-1

% Next we populate the off-diagonals, i.e. a_E and a_W
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

temperatures = horzcat(T(1),transpose(temp(:)),T(end));

fprintf('Using a non-uniform grid:\n')
fprintf('T(0) = %.1f [K] at x=%.3f\n',T(1),nodes(1))
for i = 1:N
   fprintf('T(%i) = %.1f [K] at x=%.3f\n',i,temp(i),nodes(i+1))
end
fprintf('T(%i) = %.1f [K] at x=%.3f\n',(N+1),T(end),nodes(end))

plot(nodes,temperatures,'-*b')



