function [t,y] = Euler_Method(t0, y0, t_final, dt, f)
%Euler Method Approximation
%   approximates the solution to an ODE of the form dy/dt = f(t,y)
%   using the approximation y_n+1 = y_n + dt*f(t,y)
% Input:
%  t0 = start time
%  y0 = initial value
%  t_final = end time
%  dt = time-step
%  f = right side of ODE
% Output:
% t = array of times
% y = array of solutions

% pre-allocate memory
n_total = ceil( (t_final-t0)/dt ) + 1;
y = zeros(1,n_total);
t = zeros(1,n_total);

% initial conditions
t(1) = t0;
y(1) = y0;

% approximate solutions for each time-step, from t0 to t_final
for n = 1:(n_total-1)
    % check for dt overshooting t_final
    if t(n) + dt > t_final
        dt = t_final - t(n); % adjust dt to step into t_final exactly,
                             % instead of overshooting
    end
    
    % calculate slope at t(n)
    y_prime = f(t(n),y(n));
    
    % calculate the next value using Euler's method
    y(n+1) = y(n) + dt*y_prime;
    
    % calculate next time
    t(n+1) = t(n) + dt;
end

end