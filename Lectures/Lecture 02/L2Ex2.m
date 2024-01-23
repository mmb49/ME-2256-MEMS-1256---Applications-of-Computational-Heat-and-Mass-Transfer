% Lecture 2 Example #2
clear all
close all
clc

% Defining our constants:
L = 1; % [m]
Gamma_e = 1; % [W/m-K]
Gamma_w = 1; % [W/m-K]
A_e = 1; % [m^2]
A_w = A_e;

% Defining our domain:
N = 5; % number of control volumes
n = N+2; % number of nodes
delta_x_w = L/N;
delta_x_e = L/N;
dx = L/N;

% Defining our variables of interest:
T = zeros(1,(N+2));
coeff = zeros(N,N);
b = zeros(N,1);

% Defining our boundary conditions:
q_b = 100; % [W/m^2]
T(end) = 300; % [K]

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
% modified a_W expression considering the flux
coeff(1,1) = a_E;
% Note for CV 2 we use the modified expression for a_E, which considers
% delta_x_e/2
coeff(end,end) = a_W + (A_e*((2*Gamma_e)/delta_x_e));

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
b(1) = A_w*q_b;
b(end) = (A_e*((2*Gamma_e)/dx))*T(end);

% Obtaining the solution through matrix inversion
temp = coeff\b;

% To solve for the temperature on the boundary of the domain, we must
% employ a basic differencing scheme
T(1) = temp(1) + (q_b/Gamma_w)*(delta_x_w/2);

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

x_analytic = linspace(0,L,100);
T_analytic = @(x) 300 + (q_b/Gamma_w)*(L - x); 

plot(nodes,[T(1);temp(:);T(end)],'-*b')
hold on
plot(x_analytic,T_analytic(x_analytic),'k')
ylabel('Temperature [K]')
xlabel('Distance [m]')

