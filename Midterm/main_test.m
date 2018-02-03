close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the main file for part 3 of the midterm assignment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initial conditions
y_start = [0; 1; 0; 0; 2; 1; ];
t_start = 0;
t_final = 1;
dt = [0.1 0.05 0.025 0.0125 0.00625];

%vector-function definition of the system of ODE
f_approx = @(t,y) [
    y(2);
    sin(1-exp(y(3)));
    1/(1 - log(y(5) - 1));
    (0.5*t^2 + t)*y(2) + exp(y(3))*y(1);
    1 - y(5);
    -3*((exp(y(3)) - 1)/(1+t^3))^2;
];

%vector-function for exact solution at time t
f_exact = @(t) [
  sin(t);
  cos(t);
  log(1+t);
  (0.5*t^2 + t)*sin(t);
  1 + exp(-t);
  1/(1+t^3);
];

%pre-allocate memory
max_error = zeros(1,length(dt));

for i= 1:length(dt)
    
	%test the function which solves a system of ODE with RK4
	[t,y_approx] = solve_system_RK4(f_approx,y_start,t_start,t_final,dt(i));

	%get exact solutions
	y_exact = solve_system_exact(f_exact, t, y_start);

	%calculate max_error
    error = zeros(1,length(t));
    for n= 1:length(t)
        sum = 0;
        for k= 1:length(y_start)
            sum = sum + (y_exact(k,n) - y_approx(k,n))^2;
        end
        error(n) = sqrt(sum);
    end
    
    %get the max error for each dt trial
    max_error(i) = max(error);
end

%calculate order of accuracy
order_k = zeros(1,length(dt));
for i= 1:(length(dt)-1)
   order_k(i+1) = log(max_error(i)/max_error(i+1))/log(dt(i)/dt(i+1));
end

%print table of errors and order of accuracy k
fprintf("\tdt\t\t&\tMax Error\t&\tOrder k \\\\\n");
for i= 1:(length(dt))
    fprintf("%f\t&\t%g\t&\t%f \\\\\n", dt(i), max_error(i), order_k(i));
end