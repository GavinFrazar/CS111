close all;
clear;
clc;

%parameters
t_start = 0;
t_final = 0.2;
C = 0.75;
xL = -1;
xR = 1;
yB = -1;
yT = 1;
r = 0.2;
R = 0.5;
epsilon = 0.1;
Nx_array = [50, 100, 200, 400];
Ny_array = [50, 100, 200, 400];
num_trials = length(Nx_array);
error = zeros(1,num_trials);

%define functions
f = @(t,x,y) 0;
vel = @(t,x,y) meshgrid(-y,x);
boundary_condition = @(t,x,y) 0;
c_start = @(x,y) 0.5*(1-tanh( (sqrt((x-R).^2 + y.^2)-r)/epsilon ));
c_exact = @(t,x,y) 0.5*(1-tanh( (sqrt((x-R*cos(t)).^2+(y-R*sin(t)).^2)-r)/epsilon ));

for trial = 1:num_trials
    Nx = Nx_array(trial);
    Ny = Ny_array(trial);
    x = linspace(xL,xR,Nx)';
    y = linspace(yB,yT,Ny);
    sol = upwind_scheme(t_start,t_final,C,x,y,Nx,Ny,f,boundary_condition,vel,c_start);
    exact_solution = c_exact(t_final, x, y);
    error(trial) = max(max(abs((exact_solution-sol))));
end

%calculate order of accuracy
order = zeros(1,num_trials);
for trial = 1:(num_trials-1)
    order(trial+1) = log(error(trial)/error(trial+1))/log(Nx_array(trial+1)/Nx_array(trial));
end

%print results
fprintf("\tN\t\t&\tMax Error\t&\tOrder k\t\t&\\\\\n");
for i = 1:num_trials
    fprintf("%f\t&\t%f\t&\t%f\t&\\\\\n",Nx_array(i),error(i),order(i));
end