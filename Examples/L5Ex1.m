% Lecture 5 Example #1

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


% Gaussian Elimination
AB = [coeff,b];
[m, n] = size(AB);

tStart = tic;

for j = 1:m-1 %
    for i = j+1:m 
        AB(i,j:n) = AB(i,j:n) - AB(i,j)/AB(j,j)*AB(j,j:n); % this does the lower triangular elimination
    end
end

DD = AB;

for j = 1:m-1
    for i = 1:m-j 
        DD(i,i+j:n) = DD(i,i+j:n) - DD(i,i+j)/DD(i+j,i+j)*DD(i+j,i+j:n);
    end
end

temp = zeros;
for i = 1:m
    temp(i) = DD(i,n)/DD(i,i);
end

tEnd = toc(tStart);

T(1) = 100; % [K]
T(end) = 200; % [K]

fprintf('Our calculated solution is as follows:\n')
fprintf('T(0) = %.1f [K]\n',T(1))
for i = 1:N
   fprintf('T(%i) = %.1f [K]\n',i,temp(i))
end
fprintf('T(%i) = %.1f [K]\n',(N+1),T(end))
fprintf('This took %.2e [s]\n\n',tEnd)

disp('Using matrix inversion:')
tStart = tic;
z = coeff\b;
tEnd = toc(tStart);

fprintf('T(0) = %.1f [K]\n',T(1))
for i = 1:N
   fprintf('T(%i) = %.1f [K]\n',i,z(i))
end
fprintf('T(%i) = %.1f [K]\n',(N+1),T(end))
fprintf('This took %.2e [s]\n\n',tEnd)



