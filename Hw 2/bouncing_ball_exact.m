function [x, y, vx, vy] = bouncing_ball_exact(t,a,b,c,d,r,g,alpha,beta,x0,y0,vx0,vy0)
% Calculates the exact position of the ball and its velocity at time t
% assuming the ball hits three boundaries of the container in the order:
% right, bottom, left
%   input:
%       t = time
%       a,b,c,d = locations of the left, right, bottom and top boundaries
%       of the container
%       r = ball radius
%       g = free-fall acceleration
%       alpha, beta = normal and tangent damping coefficients
%       x0, y0, vx0, vy0 = initial conditions

t0 = 0;

t1 = (b-r-x0)/vx0;

if (t < t1)
    x = x0 + vx0*(t - t0);
    y = y0 + vy0*(t - t0) - g*(t-t0)^2/2;
    vx = vx0;
    vy = vy0-g*(t-t0);
    
    return;
end

x1 = b-r;
y1 = y0 +vy0*t1-g*t1^2/2;
vx1 = -alpha*vx0;
vy1 = (vy0-g*t1)*beta;

t2 = t1 + (vy1+sqrt(vy1^2+2*g*(y1-r-c)))/g;

if (t < t2)
    x = x1 + vx1*(t - t1);
    y = y1 + vy1*(t - t1) - g*(t-t1)^2/2;
    vx = vx1;
    vy = vy1-g*(t-t1);
    
    return;
end

x2 = x1+vx1*(t2-t1);
y2 = r+c;
vx2 = beta*vx1;
vy2 = -alpha*(vy1-g*(t2-t1));

t3 = t2 + (a+r-x2)/vx2;

if (t < t3)
    x = x2 + vx2*(t - t2);
    y = y2 + vy2*(t - t2) - g*(t-t2)^2/2;
    vx = vx2;
    vy = vy2-g*(t-t2);
    
    return;
end

x3 = a+r;
y3 = y2+vy2*(t3-t2)-g*(t3-t2)^2/2;
vx3 = -alpha*vx2;
vy3 = beta*(vy2-g*(t3-t2));

x = x3 + vx3*(t - t3);
y = y3 + vy3*(t - t3) - g*(t-t3)^2/2;
vx = vx3;
vy = vy3-g*(t-t3);

end