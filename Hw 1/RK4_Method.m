function [t,y] = RK4_Method(t0, y0, t_final, dt, f)
%Classic Runge-Kutta method of approximation.
%Approximates, via Runge-Kutta, the solution for an ODE of the form y' = f(t,y)
% Input:
%  t0 = start time
%  y0 = initial value
%  t_final = end time
%  dt = time-step
%  f = right side of ODE
% Output:
% t = array of times
% y = array of solutions
n_total = ceil( (t_final - t0)/dt ) + 1;
t = zeros(1, n_total);
y = zeros(1, n_total);
t(1) = t0;
y(1) = y0;

for n = 1:(n_total-1)
   % check for dt overshooting t_final
   if t(n) + dt > t_final
       dt = t_final - t(n); % adjust dt to step into t_final exactly,
                            % instead of overshooting
   end
   
   % sample four points between y_n and y_n+1
   k1 = f(t(n),y(n));
   k2 = f(t(n) + dt/2, y(n) + k1*dt/2);
   k3 = f(t(n) + dt/2, y(n) + k2*dt/2);
   k4 = f(t(n) + dt, y(n) + k3*dt);
   
   % calculate the next value
   y(n+1) = y(n) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
   
   % increase the time
   t(n+1) = t(n)+dt;
end
end