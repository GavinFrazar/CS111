function [t, y] = solve_system_RK4(f, y0, t_start, t_final, dt)
%Solves a general system of first-order ODEs of size m
%using fourth-order accurate Runge-Kutta Method.
%INPUT: f = vectored valued function,
%       y0 = vector of initial conditions
%       t_start = starting time
%       t_final = final time
%       dt = time-step

%calculate steps
n_total = ceil((t_final - t_start)/dt);

%determine size of the system
m = length(y0);

%pre-allocate memory
t = zeros(1, n_total);
y = zeros(m, n_total);

%impose initial conditions
y(:,1) = y0;
t(1) = t_start;
n = 1;

tolerance = exp(-16);
while ( abs(t(n) - t_final) > tolerance )
   % check for dt overshooting t_final
   if t(n) + dt > t_final
       dt = t_final - t(n); % adjust dt to step into t_final exactly,
                            % instead of overshooting
   end
   
   % sample four points between y_n and y_n+1
   k1(:,1) = f(t(n),y(:,n));
   k2(:,1) = f(t(n) + dt/2, y(:,n) + k1(:,1)*dt/2);
   k3(:,1) = f(t(n) + dt/2, y(:,n) + k2(:,1)*dt/2);
   k4(:,1) = f(t(n) + dt, y(:,n) + k3(:,1)*dt);
   
   % calculate the next value
   y(:,n+1) = y(:,n) + dt*(k1(:,1) + 2*k2(:,1) + 2*k3(:,1) + k4(:,1))/6;
   
   % increase the time
   t(n+1) = t(n)+dt;
   n = n+1;
end
end

