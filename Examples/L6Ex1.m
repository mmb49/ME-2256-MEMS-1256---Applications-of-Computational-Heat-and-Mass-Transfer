%% =======================================================================
% Coded by Bailey Cassler
% =======================================================================

clc
clear all
close all

%% =======================================================================
% Problem Statement
% =======================================================================
% Consider a 2D, transient diffusion problem. A 40 [cm] by 40 [cm] plate,
% with a thermal diffusity of 0.97 [cm/s], is at an initial temperature of
% 100 deg-C. At time t=0, the sides at x(0) and x(L) are set to 0 deg-C,
% while the lateral sides are insulated. You are to plot T(x,y,t) for the
% times given by the t vector on line 18, and determine the time it takes
% to reach a maximum temperature of 10 deg-C.


%% =======================================================================
% Analytic Solution
% =======================================================================

% Analytic solution:
x_analytic = linspace(0,40,50);                      % x-domain, in cm
y_analytic = linspace(0,40,50);                      % y-domain, in cm
alpha = 0.97;                               % alpha value 9.7e-5 [m^2/s] = 0.97 [cm/s]
%t_analytic = linspace(0,300,15);
t_analytic = [10 15 30 50 100 239.8951];             % time vector
eigen = 500;                                % number of eigenvalues
T_analytic = zeros(size(x_analytic,2),size(y_analytic,2),size(t_analytic,2)+1); % predefine temperature matrix
T_analytic(:,:,1) = 100;                             % initial condition

for j=1:size(t_analytic,2) 
    T_analytic(1,:,j+1) = 0;                         % left hand side boundary condition x(0)=0
    T_analytic(end,:,j+1) = 0;                       %right hand side boudnary condition x(0)=L
    for ix = 1:length(x_analytic)
        for iy = 1:length(y_analytic)
            for n=1:2:eigen
                for m=1:2:eigen 
                    A(m,n) = 0.25*((40*(1-cos(m*pi))/(pi*m))*(40*(1-cos(n*pi))/(pi*n)));
                    lam(m,n)=((alpha*pi)/40)*sqrt(m^2+n^2);
                    T_analytic(ix,iy,j+1)=T_analytic(ix,iy,j+1)+A(m,n)*sin(m*pi*x_analytic(ix)/40)*sin(n*pi*y_analytic(iy)/40)*exp(-lam(m,n)^2*t_analytic(j));
                end
            end
        end
    end
end

% =======================================================================
% Surface Plot
% =======================================================================
[s,p] = meshgrid(x_analytic,y_analytic);
figure(1)
for i = 1:size(T_analytic,3)
    surf(s,p,T_analytic(:,:,i)); shading interp
    h = colorbar;
    h.Label.String = 'T [^{\circ}C]';
    colormap(jet)
    xlabel('x [cm]')
    ylabel('y [cm]')
    zlabel('T [^{\circ}C]')
    title('Analytic')
    zlim([0 100])
    pause(500/1000)
end



%% =======================================================================
% Finite Volume Method
% Boundary Conditions: Zero Temperature on x-edges, Zero Temperature on y-edges
% =======================================================================

clear all


% =======================================================================
% General Parameters
% =======================================================================

% Number of CVs (in each direction)
N = 100;
N2 = N+2; %(number of points in each direction)

% Domain
x = 40; %[cm] (length in j-direction)
y = 40; %[cm] (length in i-direction)
delta_x = x/N;
delta_y = y/N;
delta_z = 1;
del_x = delta_x;
del_y = delta_y;
del_z = delta_z;

% Cell Centers
cc_x = 0+(del_x/2):delta_x:x;
cc_y = 0+(del_y/2):delta_y:y;
% Grid
for k = 2:N2
    if k > 1 && k < N2
        domain_x(k) = cc_x(k-1);
        domain_y(k) = cc_y(k-1);
    elseif k == N2
        domain_x(k) = x;
        domain_y(k) = y;
    end
end
domain_y = domain_y';
[s,p] = meshgrid(domain_x,domain_y);


% Initial Conditions
T = 100*ones(N2,N2,1);
T_x0 = 0; %T(x=0) @ t=0
T_xL = 0; %T(x=L) @ t = 0;
T(:,1,1) = T_x0;
T(:,N2,1) = T_xL;
T_y0 = 0; %T(y=0) @ t=0
T_yL = 0; %T(y=L) @ t = 0;
T(1,:,1) = T_y0;
T(N2,:,1) = T_yL;

% Diffusivity
alpha = 0.97;

% Total time
t = 300;
% Delta t for each timestep
delta_t = 0.01;
% Total timesteps
total_timesteps = t/delta_t;

% Constants
k1 = alpha*delta_t/(delta_x)^2;
k2 = alpha*delta_t/(delta_y)^2;



% =======================================================================
% Explicit Scheme
% =======================================================================
% Forward Time and Central Space (FTCS) Scheme

% Tolerance
tol = 10^-6;
% Timestep
timestep = 1;
count = 1;
% Error
max_error = 1;
abs_error = zeros(N,N);


% Loop over total number of timesteps
for k = 1:total_timesteps
    % Iterate spatially until maximum absolute error is less than tolerance
    while(max_error > tol)
        % Explicit Scheme
        for j = 2:N+1
            for i = 2:N+1
                T(i,j,timestep+1) = (1-2*k1-2*k2)*T(i,j,timestep) + k1*(T(i-1,j,timestep) + T(i+1,j,timestep)) + k2*(T(i,j-1,timestep) + T(i,j+1,timestep));
            end
        end
        % Calculate absolute error between timesteps
        for j = 1:N2
            for i = 1:N2
                abs_error(i,j) = abs(T(i,j,timestep)-T(i,j,timestep+1));
            end
        end
        % Calculate max error
        max_error = max(max(abs_error));
        % Calculate max temperature
        max_T = max(max(T(:,:,timestep)))
            % Save timestep at which max temperature is 10C
            if max_T < 10
               timestep_to_10C(count) = timestep+1;
               break;
            end
        % Increment timestep
        timestep = timestep+1;
        count = count+1;
    end
end
% Timestep to reach 10C
timestep_to_10C = min(nonzeros(timestep_to_10C));
% Total time to reach 10C
time_to_10C = timestep_to_10C*delta_t;
fprintf(['The time to reach a maximum temperature of 10C is: ', num2str(time_to_10C), ' seconds.'])



% =======================================================================
% Surface Plot
% =======================================================================
figure(2)

for i = 1:50:timestep_to_10C
    surf(s,p,T(:,:,i)); shading interp
    h = colorbar;
    h.Label.String = 'T [^{\circ}C]';
    colormap(jet)
    xlabel('x [cm]')
    ylabel('y [cm]')
    zlabel('T [^{\circ}C]')
    title('Numeric')
    zlim([0 100])
    pause(1/1000)
end

