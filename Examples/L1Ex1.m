clc;
clear all;
close;
% Set maximum number of partial sums
nmax=50;
% Initialize domain
x=linspace(0,40,1000);
t=[5 10 15 20 40 80 160 240];

T = zeros(size(x,2),size(t,2)+1);

% Initial condition
for i=1:length(x)
    if x(i)<=x(length(x)) && x(i)>0
        T(i,1)=x(i);
    else
        T(i,1)=0;
    end
end

% Calculate temperature of rod at various times

for k=1:8
    
    % set boundary conditions to zero
    T(1,k) = 0;
    T(end,k) = 0;
    
    for i = 2:length(x)-1
        for n=1:2:nmax
            T(i,k+1)=T(i,k+1)+(80/pi)*(1/n)*exp((-(n^2)*(pi^2)*t(k))/1600)*sin(n*pi*x(i)/40);
        end
        for n=2:2:nmax
            T(i,k+1)=T(i,k+1)-(80/pi)*(1/n)*exp((-n^2*pi^2*t(k))/1600)*sin(n*pi*x(i)/40);
        end
    end
    
end

% contruct thses to match the temperature matrix T, [1000x9]
% [x,t]=meshgrid(1:1:40,1:10:240);

figure(1)
for i = 1:9
    plot(x,T(:,i))
    hold on
end
legend('0 s','5 s','10 s','15 s','20 s','40 s','80 s','160 s','240 s')

% figure (3)
% surf(x,t,T)
% colorbar
% xlabel('Length along rod,x')
% ylabel('Time,sec')
% zlabel ('Temperature,C')