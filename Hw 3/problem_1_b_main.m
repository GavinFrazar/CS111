close all;
clear;
clc;

%parameters
t_start = 0;
t_final = 15;
rc = 0.2;
xc = 0;
yc = 0;
U = 0.1;
C = 0.75;
xL = -1;
xR = 1;
yB = -1;
yT = 1;
r = 0.2;
R = 0.5;
Nx = 100;
Ny = 100;
epsilon = 0.05;
x0 = -0.7;
y0 = -0.2;
alpha = [0, 0.1, 0.2];
snapshot = linspace(0,15,50);
%snapshot = [0 5 10 15]; %uncomment for problem specific snapshots
num_trials = length(alpha);
error = zeros(1,num_trials);

%define functions
f = @(t,x,y) 0;
boundary_condition = @(t,x,y) 0;
c_start = @(x,y) 0.5*(1-tanh((sqrt((x-x0).^2 + (y-y0).^2)-r)/epsilon));

for trial = 1:num_trials
    x = linspace(xL,xR,Nx)';
    y = linspace(yB,yT,Ny);
    vel = @(t,x,y) cylinder_vel_field(t,x,y,alpha(trial));
    t0 = t_start;
    c0 = c_start;
    for k = 1:length(snapshot)
        sol = upwind_scheme(t0,snapshot(k),C,x,y,Nx,Ny,f,boundary_condition,vel,c0);
        
        % Do not uncomment this unless length(snapshot) is small. Opens a
        % separate figure for each snapshot.
        %%% =========================================================================================================|
        %figure('NumberTitle', 'off', 'Name', ["alpha = " + num2str(alpha(trial)) + ", t = " + num2str(snapshot(k))]);
        %%% =========================================================================================================|
        contourf(x,y,sol','LevelList', linspace(0,1,50), 'LineColor', 'none');
        colorbar;
        hold on
        % % plot velocity arrows
        % % -----------------------------------------------------------------------
        nx_quiver = 20; % coarse grid for plotting velocity field
        ny_quiver = 20;

        x_quiver = linspace(xL, xR, nx_quiver);
        y_quiver = linspace(yB, yT, ny_quiver);

        for i = 1:nx_quiver
            for j = 1:ny_quiver
                [vx_quiver(i,j),vy_quiver(i,j)] = vel(0,x_quiver(i),y_quiver(j));
            end
        end
        quiver(x_quiver, y_quiver, vx_quiver', vy_quiver','Color','w');
        % % -----------------------------------------------------------------------
        draw_disk(xc,yc,r);
        hold off
        drawnow;
        axis([xL xR yB yT]);
        axis equal;
        t0 = snapshot(k);
        c0 = @(x,y) sol;
        fprintf("alpha = %.1f | t = %u\n", alpha(trial),snapshot(k));
    end
end

function [vx,vy] = cylinder_vel_field(t,x,y,alpha)
rc = 0.2;
U = 0.1;
    idx = (sqrt(x.^2 + y.^2) >= rc);    %logical index matrix for imposing
                                        %the velocity conditional
    vx =  idx.*( U*(1 + (rc^2)*((y.^2 - x.^2)./(x.^2 + y.^2).^2)) + 2*alpha*U*rc*(y./((x.^2 + y.^2).^(3/2))) );
    vy = idx.*( -2*U*(rc^2)*((x.*y)./((x.^2 + y.^2).^2)) - 2*alpha*U*rc*(x./((x.^2 + y.^2).^(3/2))) );
end