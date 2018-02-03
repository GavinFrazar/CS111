close all;
clear;
clc;

%parameters
xL = 0;
xR = 12;
yB = 0;
yT = 3;
lambda = 0.2;
c_limit = 0.006;    

%velocity
vx = -0.8;
vy = -0.4;

%functions
c_start =@(t,x,y) y.*x*0;
c_bc = @(t,x,y) y.*x*0;
g =@(t,x,y) y.*x*0;

%time
t_start = 0;
t_final = 10;
dt = 0.1;
n_total = ceil((t_final - t_start)/dt);
t = zeros(1,n_total);

%grid space
Nx = 160;
Ny = 40;
x = linspace(xL,xR,Nx);
dx = x(2) - x(1);
y = linspace(yB,yT,Ny).';
dy = y(2) - y(1);

%logical matrix masks
BOUNDARY_MASK = get_border(Nx,Ny);
Y_BOTTOM_MASK = zeros(Ny,Nx);
Y_BOTTOM_MASK(1,2:(Nx-1)) = 1;
BOUNDARY_MASK = xor(BOUNDARY_MASK,Y_BOTTOM_MASK);

%row and column of zeros for padding
row_pad = zeros(1,Nx);
col_pad = zeros(Ny,1);

%preallocate concentration arrays
beach_a = floor(4/dx + 1);
beach_b = floor(6/dx + 1);
beach_c = floor(8/dx + 1);
conc_a = zeros(1,n_total);
conc_b = zeros(1,n_total);
conc_c = zeros(1,n_total);

%tolerance for equality checking of floating point variables
tol = dt/10;

%impose initial conditions
A = get_A(Nx,Ny,dt,dx,dy,lambda,vy);
c_old = c_start(t_start,x,y);
c_new = c_old;
conc_a(1) = c_old(1,beach_a);
conc_b(1) = c_old(1,beach_b);
conc_c(1) = c_old(1,beach_c);
n = 1;
t(1) = t_start;
while t(n) < t_final
   if ( abs(t(n) - 1) <= tol || abs(t(n) - 4) <= tol || abs(t(n) - 7) <= tol)
       figure('NumberTitle', 'off', 'Name', "t = " + num2str(t(n)));
       contourf(x,y,c_new, 100, 'LineColor','none');
        colorbar;
        caxis([0,0.016])
        drawnow;
        axis([xL xR yB yT]);
        axis equal; 
   end
    
   if t(n) + dt > t_final
       dt = t_final - t(n);
        A = get_A(Nx,Ny,dt,dx,dy,lambda,vy);
   end
   
   t(n+1) = t(n) + dt;
   dcdx = diff([c_old, col_pad],1,2)/dx;
   dcdy = diff([c_old; row_pad],1,1)/dy;
   rhs = BOUNDARY_MASK.*(c_bc(t(n+1),x,y))... %rhs for bc
        + Y_BOTTOM_MASK.*( c_old + dt*f(t(n+1),x,y)... %rhs for robinson bc
            - vx*dt*dcdx...
            - vy*dt*dcdy...
            - (2*dt/dy)*g(t(n+1),x,y) )...
        + ~(BOUNDARY_MASK | Y_BOTTOM_MASK).*( c_old + dt*f(t(n+1),x,y)... %rhs for internal points
            - vx*dt*dcdx...
            - vy*dt*dcdy );
   rhs = reshape(rhs.',[],1); %make rhs into a column vector
   c_new = A\rhs;
   c_new = reshape(c_new,[],Ny)'; %reshape solution as a matrix
   c_old = c_new;
   
   %observe oil concentration at each beach
   conc_a(n+1) = c_new(1,beach_a);
   conc_b(n+1) = c_new(1,beach_b);
   conc_c(n+1) = c_new(1,beach_c);
   n = n+1;
end

figure('NumberTitle', 'off', 'Name', "Beach A (2,0) = ");
plot(t, conc_a, '-', 'Linewidth', 2);
hold on
close_beach = refline(0,c_limit);
close_beach.Color = 'r';
hold off
set(gca,'FontSize',18)
xlabel('Time');
ylabel('Concentration');
h = legend(close_beach, "Oil conc. = 0.006");
rect = [0.4, 0.8, .08, .1];
set(h, 'Position', rect)

figure('NumberTitle', 'off', 'Name', "Beach B (4,0) = ");
plot(t, conc_b, '-', 'Linewidth', 2);
hold on
close_beach = refline(0,c_limit);
close_beach.Color = 'r';
hold off
set(gca,'FontSize',18)
xlabel('Time');
ylabel('Concentration');
h = legend(close_beach, "Oil concentration = 0.006");
set(h,'Location','best');

figure('NumberTitle', 'off', 'Name', "Beach C (6,0) = ");
plot(t, conc_c, '-', 'Linewidth', 2);
hold on
close_beach = refline(0,c_limit);
close_beach.Color = 'r';
hold off
set(gca,'FontSize',18)
xlabel('Time');
ylabel('Concentration');
h = legend(close_beach, "Oil concentration = 0.006");
set(h,'Location','best');

function out = f(t,x,y)
epsilon = 0.1;
x_s = 10;
r_s = 0.1;
    if t > 0.5
        out = 0;
    else
        out = 0.5*(1 - tanh(( sqrt(y.^2 + (x - x_s).^2) - r_s)/epsilon));
    end
end