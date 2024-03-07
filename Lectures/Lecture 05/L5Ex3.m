% Lecture 5 Example #3

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
N = 5000; % number of control volumes
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

% Timing the array population and solution calcuation
tStart = tic;

% Populating our coefficient matrix.
coeff = diag(a_P*ones(N,1)) + diag(-a_E*ones(N-1,1),1) + diag(-a_W*ones(N-1,1),-1);

% Next we have to populate a_P for CVs 1 and N. Note for CV 1 we use the
% modified a_W expression, which considers delta_x_w/2
coeff(1,1) = a_E + (A_w*((2*Gamma_w)/dx));
% Note for CV 2 we use the modified expression for a_E, which considers
% delta_x_e/2
coeff(end,end) = a_W + (A_e*((2*Gamma_e)/dx));

% Next we populate the RHS matrix
b(1) = (A_w*((2*Gamma_w)/dx))*T(1);
b(end) = (A_e*((2*Gamma_e)/dx))*T(end);

temperature = coeff\b;
tEnd_matrix = toc(tStart);

fprintf('Using matrix inversion:\n')
%fprintf('T(0) = %.1f [K]\n',T(1))
%for i = 1:N
%   fprintf('T(%i) = %.1f [K]\n',i,temperature(i))
%end
%fprintf('T(%i) = %.1f [K]\n',(N+1),T(end))
fprintf('This took %.2e [s]\n\n',tEnd_matrix)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timing the coefficient population and solution
tStart = tic;

% TDMA - we first need to extract the diagonal and off-diagonal values
b = diag(coeff); % diagonal
a = diag(coeff,-1); % lower diagonal
c = diag(coeff,1); % upper diagonal

% redefining the RHS array to match our convection
d = zeros(1,N);
d(1) = (A_w*((2*Gamma_w)/dx))*T(1);
d(end) = (A_e*((2*Gamma_e)/dx))*T(end);

% initializing our vectors
m = zeros(N,1);
d_prime = zeros(N,1);

% Our off diagonals need to be the same size as the diagonal, thus we add
% zero to the proper index
c = [c; 0];
a = [0; a];

% We need to initialize our first m value, m(1) since our loop starts at 2
m(1) = c(1)/b(1); 
d_prime(1) = d(1)/b(1);

for i=2:N
    m(i) = c(i)/(b(i)-a(i)*m(i-1)) ;
    d_prime(i) = (d(i)-a(i)*d_prime(i-1))/(b(i)-a(i)*m(i-1)) ;
end

phi = zeros(N, 1) ;
phi(N) = d_prime(N) ;
for i=N-1:-1:1
    phi(i) = d_prime(i) - m(i)*phi(i+1);
end

tEnd_TDMA = toc(tStart);

fprintf('Using TDMA:\n')
%fprintf('T(0) = %.1f [K]\n',T(1))
%for i = 1:N
%   fprintf('T(%i) = %.1f [K]\n',i,phi(i))
%end
%fprintf('T(%i) = %.1f [K]\n',(N+1),T(end))
fprintf('This took %.2e [s]\n\n',tEnd_TDMA)

dummy = tEnd_matrix/tEnd_TDMA;

fprintf('The TDMA algorithm is %.2f times faster than matrix inversion\n',dummy)