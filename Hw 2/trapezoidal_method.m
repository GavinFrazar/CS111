function [x_error, y_error, vx_error, vy_error, x, y] = trapezoidal_method(t0, t_final, x0, y0, vx0, vy0, dt, f_x, f_y, f_vx, f_vy, a, b, c, d, r, alpha, beta, g)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n_total = ceil((t_final - t0)/dt);
n = 1;
x = zeros(1, n_total);
y = zeros(1, n_total);
vx = zeros(1, n_total);
vy = zeros(1, n_total);
out_of_bounds_x = false;
out_of_bounds_y = false;

% impose initial conditions
t = t0;
x(1) = x0;
y(1) = y0;
vx(1) = vx0;
vy(1) = vy0;

while (t < t_final)
    % check for overshooting
    if (t + dt > t_final)
        dt = abs(t_final - t);
    end
    
    t_next = t + dt;
    
    %calculate the next position and velocity vectors with the trapezoidal
    %approximation
    x_temp = x(n) + dt*f_x(t, x(n), vx(n));
    vx_temp = vx(n) + dt*f_vx(t, x(n), vx(n));
    x(n+1) = x(n) + 0.5*dt*(f_x(t, x(n), vx(n)) + f_x(t_next, x_temp, vx(n)));
    vx(n+1) = vx(n) + 0.5*dt*(f_vx(t, x(n), vx(n)) + f_vx(t_next, x_temp, vx_temp));
    
    y_temp = y(n) + dt*f_y(t, y(n), vy(n));
    vy_temp = vy(n) + dt*f_vy(t, y(n), vy(n));
    y(n+1) = y(n) + 0.5*dt*(f_y(t, y(n), vy(n)) + f_y(t_next, y_temp, vy_temp));
    vy(n+1) = vy(n) + 0.5*dt*(f_vy(t, y(n), vy(n)) + f_vy(t_next, y_temp, vy_temp));
    
    %apply consequences of collision with a wall and restore dt to its
    %value before the collision.
    if (out_of_bounds_x)
       dt = dt_prev;
       out_of_bounds_x = false;
       vx(n+1) = -vx(n+1)*alpha;
       vy(n+1) = vy(n+1)*beta;
    end
    
    if (out_of_bounds_y)
       dt = dt_prev;
       out_of_bounds_y = false;
       vy(n+1) = -vy(n+1)*alpha;
       vx(n+1) = vx(n+1)*beta;
    end
    
    % check for bounds
    %left wall
    if (x(n+1) < a+r)
        dt_prev = dt;
        dt = dt*(abs(a+r-x(n)))/abs(x(n+1)-x(n));
        out_of_bounds_x = true;
        t = t+dt-dt_prev;
        continue;
    end
    %right wall
    if (x(n+1) > (b-r))
        dt_prev = dt;
        dt = dt*(abs(b-r-x(n)))/abs(x(n+1)-x(n));
        out_of_bounds_x = true;
        t = t+dt-dt_prev;
        continue;
    end
    %floor
    if (y(n+1) < (c+r))
        dt_prev = dt;
        dt = dt*(abs(c+r-y(n)))/abs(y(n+1)-y(n));
        out_of_bounds_y = true;
        t = t+dt - dt_prev;
        continue;
    end
    %ceiling
    if (y(n+1) > (d-r))
        dt_prev = dt;
        dt = dt*(abs(d-r-y(n)))/abs(y(n+1)-y(n));
        out_of_bounds_y = true;
        t = t+dt-dt_prev;
        continue;
    end
    
    %calculate exact values at this time
    [x_exact, y_exact, vx_exact, vy_exact] = bouncing_ball_exact(t_next,a,b,c,d,r,g,alpha,beta,x0,y0,vx0,vy0);
    %
        
    %get the errors to return
    if(t_next >= t_final)
       x_error = abs(x(n+1)-x_exact);
       y_error = abs(y(n+1)-y_exact);
       vx_error = abs(vx(n+1)-vx_exact); 
       vy_error = abs(vy(n+1)-vy_exact);
       
       %{
       fprintf("\nTime: %f", t_next);
       fprintf("\nx_exact = %f, x_approx = %f", x(n+1), x_exact);
       fprintf("\ny_exact = %f, y_approx = %f", y(n+1), y_exact);
       fprintf("\nvx_exact = %f, vx_approx = %f", vx(n+1), vx_exact);
       fprintf("\nvy_exact = %f, vy_approx = %f\n", vy(n+1), vy_exact);
       %}
    end
    
    
    %increase time
    t = t + dt;
    n = n + 1;
end

end

