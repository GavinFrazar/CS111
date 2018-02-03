close all;
clear;
clc;

%parameters
%space (meters)
xL = -2;
xR = 2;
yB = -2.5;
yT = 2.5;

%grid resolution
Nx = 80;
Ny = 100;

%time (seconds)
t_start = 0;
t_final = 1500;
dt = 5;

%thermal diffusivity
lambda = 1.5 * 10^(-3);

%temperature (celsius)
%boundary conditions
T_bc = @(t,x,y) ones(Ny,Nx).*min(20+4*t/3,100);

%initial conditions
T_start =@(t,x,y) 20;

%source term
f = @(t,x,y) 0;

%cooking temp
T_cooking = 65;

%discretization in space
x = linspace(xL,xR,Nx);
y = linspace(yB,yT,Ny)';
dx = x(2) - x(1);
dy = y(2) - y(1);

%preallocate
temp = zeros(1,ceil((t_final - t_start)/dt));
t = zeros(1, ceil((t_final - t_start)/dt));
center_x = ceil(Nx/2);
center_y = ceil(Ny/2);

%logical matrix mask
BOUNDARY = get_border(Nx,Ny);

%variables to keep track of temperature of the center of the potato
t_began_cooking = 0;
already_cooking = false;
temp(1) = T_start(t_start,x(center_x),y(center_y));

%impose initial conditions
A = get_A(Nx,Ny,dt,dx,dy,lambda);
T_old = ones(Ny,Nx).*T_start(t_start,x,y);
T_new = T_old;

t(1) = t_start;
n = 1;
while t(n) < t_final
    if (t(n) == 0 || t(n) == 200 || t(n) == 400 || t(n) == 600)
        %thermal image
        figure('NumberTitle', 'off', 'Name', "t = " + num2str(t(n)));
        contourf(x,y,T_new, 100, 'LineColor','none');
        colorbar;
        caxis([20,100])
        axis([xL xR yB yT]);
        axis equal;
    end

    if (temp(n) >= T_cooking) && (already_cooking == false)
        t_began_cooking = t(n);
        already_cooking = true;
    end
    
    if t(n) + dt > t_final
        dt = t_final - t(n);
        A = get_A(Nx,Ny,dt,dx,dy,lambda);
    end

    t(n+1) = t(n) + dt;
    rhs = BOUNDARY.*(T_bc(t(n+1),x,y)) + ~BOUNDARY.*(T_old + dt*f(t(n+1),x,y));
    rhs = reshape(rhs.',[],1);
    T_new = A\rhs;
    T_new = reshape(T_new,Nx,Ny)';
    temp(n+1) = T_new(ceil(Ny/2),ceil(Nx/2));
    T_old = T_new;
    n = n+1;
end

figure;
plot(t, temp, '-', 'Linewidth', 2);
hold on
cooking_threshold = refline(0,T_cooking);
cooking_threshold.Color = 'r';
hold off
set(gca,'FontSize',18)
xlabel('Time');
ylabel('Temperature');
legend(cooking_threshold, "Cooking Temperature = " + num2str(T_cooking) + " C");
fprintf("The potato started cooking at %d seconds, and finished cooking 300 seconds later at %d seconds."...
    ,t_began_cooking, t_began_cooking+300);