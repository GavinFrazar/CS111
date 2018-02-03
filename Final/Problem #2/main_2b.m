close all;
clear;
clc;

%parameters
xL = -1;
xR = 3;
yB = -1.5;
yT = 1.5;
lambda = 0.7;

%velocity
vx = -0.8;
vy = -0.4;

%functions
c_exact =@(t,x,y) exp(-t).*cos(y).*sin(x);
c_start =@(t,x,y) c_exact(0,x,y);
c_bc = @(t,x,y) c_exact(t,x,y);
f =@(t,x,y) exp(-t)*(   -cos(y).*sin(x)...
                        +vx*cos(y).*cos(x)...
                        -vy*sin(y).*sin(x)...
                        +2*lambda*cos(y).*sin(x)...
                    );
g =@(t,x,y) exp(-t)*(  -lambda*sin(y).*sin(x)...
                        -vy*cos(y).*sin(x)...
                     );

%time
t_start = 0;
t_final = 1;

%grid space
Nx_array = [20 40 80 160];
Ny_array = [15 30 60 120];

num_trials = length(Nx_array);
for trial = 1:num_trials
    %discretization of space
    Nx = Nx_array(trial);
    Ny = Ny_array(trial);
    x = linspace(xL,xR,Nx);
    dx = x(2) - x(1);
    y = linspace(yB,yT,Ny).';
    dy = y(2) - y(1);
    dt = dx/2;
    
    %logical matrix masks
    Y_BOTTOM_MASK = zeros(Ny,Nx);
    Y_BOTTOM_MASK(1,2:(Nx-1)) = 1;
    BOUNDARY_MASK = get_border(Nx,Ny);
    BOUNDARY_MASK = xor(BOUNDARY_MASK,Y_BOTTOM_MASK);
    
    %row and column of zeros to pad with
    row_pad = zeros(1,Nx);
    col_pad = zeros(Ny,1);
    
    %impose initial conditions
    A = get_A(Nx,Ny,dt,dx,dy,lambda,vy);
    c_old = c_start(t_start,x,y);
    c_new = c_old;
    t = t_start;
    while t < t_final
       if t + dt > t_final
          dt = t_final - t;
          A = get_A(Nx,Ny,dt,dx,dy,lambda,vy);
       end
       
       %
       dcdx = diff([c_old, col_pad],1,2)/dx;
       dcdy = diff([c_old; row_pad],1,1)/dy;
       rhs = BOUNDARY_MASK.*( c_bc(t+dt,x,y) )... %rhs for bc
            + Y_BOTTOM_MASK.*( c_old + dt*f(t+dt,x,y)... %rhs for robinson bc
                - vx*dt*dcdx...
                - vy*dt*dcdy...
                - (2*dt/dy)*g(t+dt,x,y) )...
            + ~(BOUNDARY_MASK | Y_BOTTOM_MASK).*( c_old + dt.*f(t+dt,x,y)... %rhs for internal points
                - vx.*dt.*dcdx...
                - vy.*dt.*dcdy );
       rhs = reshape(rhs.',[],1); %make rhs into a column vector
       %solve matrix equation for c_new
       c_new = A\rhs;
       c_new = reshape(c_new,Nx,Ny)'; %reshape solution as a matrix
       
       %increment time and step
       c_old = c_new;
       t = t+dt;
    end
   
   %compare with exact solution at t = t_final
   exact_sol = c_exact(t_final,x,y);
   error(trial) = max(max(abs(c_new - exact_sol)));
end

%calculate order of accuracy
Order_k = zeros(1,num_trials);
for trial = 1:(num_trials -1)
	Order_k(trial+1) = log(error(trial)/error(trial+1))/abs(log(Nx_array(trial+1)/Nx_array(trial)));
end

%print results
fprintf("(Nx, Ny)\t&\tError\t\t&\tOrder k\t\t\\\\\n");
for i = 1:num_trials
	fprintf("(%d, %d)\t&\t%f\t&\t%f\t \\\\\n",Nx_array(i),Ny_array(i),error(i),Order_k(i));
end