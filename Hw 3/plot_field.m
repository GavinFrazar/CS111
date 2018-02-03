% Example demonstrating plotting of 2D fields and velocity fields

clear;

xL =-1;
xR = 1;
yB =-1;
yT = 1;

% initial conditions
xi = -.7;
yi = -.2;
ri = 0.2;
eps = 0.05;

c0 = @(x,y) 0.5*(1-tanh((sqrt((x-xi)^2+(y-yi)^2)-ri)/eps));

% velocity field
vx_func = @(t,x,y) -y;
vy_func = @(t,x,y) x;

% space discretization
nx = 100;
ny = 100;

x = linspace(xL, xR, nx);
y = linspace(yB, yT, ny);


% initial conditions
for i = 1:nx
    for j = 1:ny
        c_old(i,j) = c0(x(i),y(j));
    end
end


contourf(x,y,c_old','LevelList', linspace(0.0,1.0,50), 'LineColor', 'none');
colorbar;
hold on
% % uncomment to plot velocity arrows
% % -----------------------------------------------------------------------
% nx_quiver = 20; % coarse grid for plotting velocity field
% ny_quiver = 20;
% 
% x_quiver = linspace(xL, xR, nx_quiver);
% y_quiver = linspace(yB, yT, ny_quiver);
% 
% for i = 1:nx_quiver
%     for j = 1:ny_quiver
%         vx_quiver(i,j) = vx_func(0,x_quiver(i),y_quiver(j));
%         vy_quiver(i,j) = vy_func(0,x_quiver(i),y_quiver(j));
%     end
% end
% quiver(x_quiver, y_quiver, vx_quiver', vy_quiver','Color','w');
% % -----------------------------------------------------------------------
hold off
axis([xL xR yB yT]);
axis equal;
