clear
clc
close all

% Use your code to approximate the derivative of the following expression at x=792.93

format long

yforward = zeros;
ybackward = zeros;
ycentral = zeros;
ye = zeros;
yf = zeros;
yb = zeros;
yc = zeros;
ef = zeros;
eb = zeros;
ec = zeros;
H = zeros;
% calculatingt the actual derivative

syms x y
y = log(sin(x)^cos(x))/x^(1/3);
dy = diff(y,'x');

disp('Using Matlab, the derivative is calculated as:')
disp(dy)

yy = (cos(792.93)^2*sin(792.93)^(cos(792.93) - 1) - log(sin(792.93))*sin(792.93)*sin(792.93)^cos(792.93))/(792.93^(1/3)*sin(792.93)^cos(792.93)) - log(sin(792.93)^cos(792.93))/(3*792.93^(4/3));
disp('The solution to the Matlab derivative |x=792.93 is:')
disp(yy)

x = 792.93;
yDexact = (((cos(x)^2*sin(x)^(cos(x) - 1)) - (log(sin(x))*sin(x)*sin(x)^cos(x))/(x^(1/3)*sin(x)^cos(x)))) - log(sin(x)^cos(x))/(3*x^(4/3));

h = 792; % this will be treated as our delta x, with an initial guess of 1
for i = 1:400
    h =  h/1.1;
    H(i) = h;
    
    x0 = 792.93; % we can treat this as x_i
    x_1 = x0 - h; % we can treat this as x_i - 1 for the backward difference
    x1 = x0 + h; % and we can treat this as x_1 + 1 for the forward
    
    y_1 = log(sin(x_1)^cos(x_1))/(x_1^(1/3)); % this is f(x_i - 1) 
    y1 = log(sin(x1)^cos(x1))/(x1^(1/3)); % this is f(x_i + 1)
    y0 = log(sin(x0)^cos(x0))/(x0^(1/3)); % this is f(x_i)
    
    yforward(i) = (y1 - y0)/h;
    ybackward(i) = (y0 - y_1)/h;
    ycentral(i) = (y1 - y_1)/(2*h);
    
    ye(i) = ((cos(x0)^2*sin(x0)^(cos(x0) - 1) - log(sin(x0))*sin(x0)*sin(x0)^cos(x0))/(x0^(1/3)*sin(x0)^cos(x0))) - log(sin(x0)^cos(x0))/(3*x0^(4/3));
    yf(i) = (y1 - y0)./h;
    yb(i) = (y0 - y_1)./h;
    yc(i) = (y1 - y_1)./(2*h);
    
    ef(i) = abs(yf(i) - ye(i));
    eb(i) = abs(yb(i) - ye(i));
    ec(i) = abs(yc(i) - ye(i));
end

disp('The forward, backward and central difference derivative values calculated via the code are:')
yfor = yforward(ef == min(ef));
ybac = ybackward(eb == min(eb));
ycen = ycentral(ec == min(ec));
disp(yfor)
disp(ybac)
disp(ycen)

figure
a = logspace(0,-12,12);
loglog(H,ef,'--k',H,eb,':g',H,ec,'-r')
axis([10E-12 10E0 10E-15 10E0])
title('\Deltax versus log(e_t)')
% set(gca,'FontSize',20)
legend('Forward Difference','Backward Difference','Central Difference')
xlabel('\Deltax')
ylabel('log(e_t)')

%% Finding the minimum delta x

d_x_min_forward = H(ef == min(ef));
d_x_min_backward = H(eb == min(eb));
d_x_min_center = H(ec == min(ec));

min_forward = [d_x_min_forward, min(ef)];
min_backward = [d_x_min_backward, min(eb)];
d_x_min_center = [d_x_min_center, min(ec)];

disp('The optimum delta_x for the forward difference method is:')
disp(d_x_min_forward)
disp('The optimum delta_x for the backward difference method is:')
disp(d_x_min_backward)
disp('The optimum delta_x for the center difference method is:')
disp(d_x_min_center)
