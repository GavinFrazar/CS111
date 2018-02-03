close all;
clear;
clc;

%spacial dimensions
xL = -1;
xR = 1;
yB = -0.5;
yT = 1.7;

%grid resolution
Nx_array = [25,50,100];
Ny_array = [30,60,120];
%time (seconds)
t_start = 0;
t_final = 1;

%thermal diffusivity
lambda = 0.75;

%temperature (celsius)
%exact solution
T_exact =@(t,x,y) cos(y).*sin(x)*exp(-t);

%boundary conditions
T_bc = @(t,x,y) T_exact(t,x,y);

%initial conditions
T_start =@(t,x,y) T_exact(t,x,y);

%source term
f = @(t,x,y) ((2*lambda - 1)*exp(-t))*cos(y).*sin(x);

num_trials = length(Nx_array);
for trial = 1:num_trials
    %discretization of space
    Nx = Nx_array(trial);
    Ny = Ny_array(trial);
    x = linspace(xL,xR,Nx_array(trial));
    y = linspace(yB,yT,Ny_array(trial))';
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    %discretization of time
    dt = dx/2;
    
    %logical matrix mask
    BOUNDARY = get_border(Nx,Ny);
    
    %impose initial conditions
    A = get_A(Nx,Ny,dt,dx,dy,lambda);
	T_old = T_start(t_start,x,y);
    T_new = T_old;
    t = t_start;
    
    while(t < t_final)
        if t + dt > t_final
            dt = t_final - t;
            
            %recalculate matrix A using new dt
            A = get_A(Nx,Ny,dt,dx,dy,lambda);
        end
        
        %create rhs of the matrix equation. Logical mask used to place
        %boundary conditions in appropriate place. Elsewhere we have the
        %value rhs = T_old + dt*f
        rhs = BOUNDARY.*(T_bc(t+dt,x,y)) + ~BOUNDARY.*(T_old + dt*f(t+dt,x,y));
        %convert rhs from matrix into column vector
        rhs = reshape(rhs.',[],1);
        
        %solve matrix equation for T_new
        T_new = A\rhs;
        
        %convert T_new from column vector back into a matrix
        T_new = reshape(T_new,Nx,Ny)';
        
        %increment time and step
        T_old = T_new;
        t = t+dt;
    end
    exact_sol = T_exact(t_final,x,y);
	error(trial) = max(max(abs(T_new - exact_sol)));
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