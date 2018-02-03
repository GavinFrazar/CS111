function [t,y] = Trapezoidal_Method(t0, y0, t_final, dt, f)
%Trapezoidal method of approximation.
%Approximates, via trapezoids, the solution for an ODE of the form y' = f(t,y)
% Input:
%  t0 = start time
%  y0 = initial value
%  t_final = end time
%  dt = time-step
%  f = right side of ODE
% Output:
% t = array of times
% y = array of solutions

%pre-allocate memory
n_total = ceil((t_final-t0)/dt) + 1;
t = zeros(1,n_total);
y = zeros(1,n_total);

% initial conditions
y(1) = y0;
t(1) = t0;

% approximate solutions for each time-step, from t0 to t_final
for n = 1:(n_total-1)
    % check for dt overshooting t_final
    if t(n) + dt > t_final
        dt = t_final - t(n); % adjust dt to step into t_final exactly,
                             % instead of overshooting
    end
    
    % Since the trapezoidal method is an implicit method of approximation,
    % we must use Euler's method at each step to temporarily calculate the
    % next value, then use that value to get calculate a (usually) better
    % approximation with the trapezoidal method.
    y_temp = y(n) + f(t(n),y(n))*dt;
    
    % plug the temp value for y_n+1 into the trapezoidal method formula
    y(n+1) = y(n) + 0.5*dt*( f(t(n),y(n)) + f(t(n+1),y_temp) );
    
    % increase the time
    t(n+1) = t(n) + dt;
end

end