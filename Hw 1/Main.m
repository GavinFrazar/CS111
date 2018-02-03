close all;
clear;
clc;

% parameters %
g = 9.81;   % Force of gravity (m/s^2)
m = 75;     % mass (kg)
c_d = 0.25; % drag coefficient (kg/m)

% initial conditions %
t0 = 0;         %intial time (seconds)
v0 = 0;         %initial velocity (m/s)
t_final = 15;   %final time (seconds)
dt = [3.3 0.1 0.05 0.025 0.0125];   %time-step trials (seconds)

% data containers
num_methods = 3; % number of approximation methods used
num_trials = length(dt);
f_parachute = @(t,v) g-(c_d/m)*v^2;   %define the parachute function
v_approx_array = cell(num_trials,num_methods);  %array of approx. velocity arrays
v_exact_array = cell(num_trials,1);   %array of exact velocity arrays
t_array = cell(num_trials,1);         %array of time arrays
errors = cell(num_trials,num_methods);%2D array of error arrays
error_maximums = zeros(num_trials, num_methods); %array of maximum error in each trial for each method

% solve the ODE using each method of approximation.
for i = 1:length(dt)
    [t_array{i,1}, v_approx_array{i,1}] = Euler_Method(t0, v0, t_final, dt(i), f_parachute); 
    [t_array{i,1}, v_approx_array{i,2}] = Trapezoidal_Method(t0, v0, t_final, dt(i), f_parachute);
    [t_array{i,1}, v_approx_array{i,3}] = RK4_Method(t0, v0, t_final, dt(i), f_parachute);
end

% solve the ODE exactly, using the analytical method, and calculate the
% error of the approximation methods at each step.
for i = 1:num_trials
   n_total = length(t_array{i,1});
   v_exact = zeros(1, n_total);
   for k = 1:num_methods
        errors{i,k} = zeros(1,n_total);
   end
   
   for j = 1:n_total
       v_exact(j) = sqrt(g*m/c_d)*tanh(t_array{i,1}(j)*sqrt(g*c_d/m));
       for k = 1:num_methods
           errors{i,k}(j) = abs(v_exact(j) - v_approx_array{i,k}(j));
       end
   end
   v_exact_array{i,1} = v_exact;
end

%Store the max error of each ith trial for each kth method of approximation
for i = 1:num_trials
    for k = 1:num_methods
        error_maximums(i,k) = max(errors{i,k}); 
    end
end

%calculate the order of accuracy when reducing dt
order = zeros(num_trials,num_methods);

for i = 2:num_trials
   for k = 1:num_methods
       order(i,k) = log(error_maximums(i-1,k)/error_maximums(i,k))/(log(dt(i-1)/dt(i)));
   end
end
fprintf('EULER METHOD\n============\n');
fprintf('\nThe max error in each trial using Euler method was: ');

%Print Euler maximum error for each time-step trial
fprintf('\n\tdt\t|\tMax_error\t|\tOrder k');
fprintf('\n%.4f\t|\t%g\t\t|\tN/A', dt(1), error_maximums(1,1)); 
for i = 2:num_trials
   fprintf('\n%.4f\t|\t%g\t|\t%g', dt(i), error_maximums(i,1), order(i,1)); 
end

fprintf('\n==================\nTRAPEZOIDAL METHOD\n==================');
fprintf('\n\nThe max error in each trial using trapezoidal method was: ');

%Print trapezoidal maximum error for each time-step trial
fprintf('\n\tdt\t|\tMax_error\t|\tOrder k');
fprintf('\n%.4f\t|\t%g\t\t|\tN/A', dt(1), error_maximums(1,2)); 
for i = 2:num_trials
   fprintf('\n%.4f\t|\t%g\t|\t%g', dt(i), error_maximums(i,2), order(i,2)); 
end

fprintf('\n==================\nRK4 METHOD\n==================');
fprintf('\n\nThe max error in each trial using RK4 was: ');

%Print RK4 maximum error for each time-step trial
fprintf('\n\tdt\t|\tMax_error\t|\tOrder k');
fprintf('\n%.4f\t|\t%g\t|\tN/A', dt(1), error_maximums(1,3)); 
for i = 2:num_trials
   fprintf('\n%.4f\t|\t%g\t|\t%g', dt(i), error_maximums(i,3), order(i,3)); 
end

%Plot results of each method for dt = 3.3, and the exact solution.
plot(t_array{1,1}, v_approx_array{1,1}, 'o', 'Linewidth', 2);
hold on;
plot(t_array{1,1}, v_approx_array{1,2}, '^', 'Linewidth', 2);
plot(t_array{1,1}, v_approx_array{1,3}, 's', 'Linewidth', 2);
plot(t_array{1,1}, v_exact_array{1,1}, '-');
hold off;
set(gca,'FontSize',18)
xlabel('time');
ylabel('velocity');
legend('Euler', 'Trapezoidal', 'RK4', 'Exact');
