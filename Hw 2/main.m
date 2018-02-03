close all;
clear;
clc;

% parameters %
g = 9.81;
r = 0.05;
a = 0;
b = 1;
c = 0;
d = 1;

% Initial conditions %
t0 = 0;
t_final = [0.2 1.0 2.5];
error_comparison_time = 0.931;
x0 = 0.1;
y0 = 0.7;
vx0 = 3;
vy0 = 1;
alpha = 0.8;
beta = 0.9;
dt_list = [0.02 0.01 0.005 0.0025 0.00125 0.000625];

% function definitions %
f_x = @(t, x, vx) vx;
f_vx = @(t, x, vx) 0;
f_y = @(t, y, vy) vy;
f_vy = @(t, y, vy) -g;

%arrays to store error of each time step at t==0.931
x_error = zeros(1,length(dt_list));
y_error = zeros(1,length(dt_list));
vx_error = zeros(1,length(dt_list));
vy_error = zeros(1,length(dt_list));

order_x = zeros(1,length(dt_list));
order_y = zeros(1,length(dt_list));
order_vx = zeros(1,length(dt_list));
order_vy = zeros(1,length(dt_list));

%solve for t_final = 0.931, using each dt
for i= 1:length(dt_list)
    [x_error(i), y_error(i), vx_error(i), vy_error(i), ~,~] = trapezoidal_method(t0, error_comparison_time, x0, y0, vx0, vy0, dt_list(i), f_x, f_y, f_vx, f_vy, a, b, c, d, r, alpha, beta, g);
end

%computer order of accurracy
for i=1:(length(dt_list)-1)
    order_x(i+1) = log(x_error(i)/x_error(i+1))/log(dt_list(i)/dt_list(i+1));
    order_y(i+1) = log(y_error(i)/y_error(i+1))/log(dt_list(i)/dt_list(i+1));
    if (vx_error(i+1) ~= 0)
        order_vx(i+1) = log(vx_error(i)/vx_error(i+1))/log(dt_list(i)/dt_list(i+1));
    end
    order_vy(i+1) = log(vy_error(i)/vy_error(i+1))/log(dt_list(i)/dt_list(i+1)); 
end

%get average of the orders of accuracy
avg_order_x = 0;
avg_order_y = 0;
avg_order_vx = 0;
avg_order_vy = 0;

for i=1:length(dt_list)
   avg_order_x = avg_order_x + order_x(i);
   avg_order_y = avg_order_y + order_y(i);
   avg_order_vx = avg_order_vx + order_vx(i);
   avg_order_vy = avg_order_vy + order_vy(i);
end
avg_order_x = avg_order_x/5;
avg_order_y = avg_order_y/5;
avg_order_vx = avg_order_vx/5;
avg_order_vy= avg_order_vy/5;

%print error and order of accurracy results
fprintf("\tdt\t\t&\tx_error\t\t&\tOrder k");
for i= 1:(length(dt_list))
    fprintf("\n%f\t&\t%g\t&\t%f", dt_list(i), x_error(i), order_x(i));
end
fprintf("\nAverage order: %f", avg_order_x);

fprintf("\n------------------------------------------");

fprintf("\n\tdt\t\t&\ty_error\t\t&\tOrder k");
for i= 1:(length(dt_list))
    fprintf("\n%f\t&\t%g\t&\t%f", dt_list(i), y_error(i), order_y(i));
end
    fprintf("\nAverage order: %f", avg_order_y);
fprintf("\n------------------------------------------");

fprintf("\n\tdt\t\t&\tvx_error\t\t&\tOrder k");
for i= 1:(length(dt_list))
    fprintf("\n%f\t&\t%g\t&\t%f", dt_list(i), vx_error(i), order_vx(i));
end
fprintf("\nAverage order: %f", avg_order_vx);
fprintf("\n------------------------------------------");

fprintf("\n\tdt\t\t&\tvy_error\t\t&\tOrder k");
for i= 1:(length(dt_list))
    fprintf("\n%f\t&\t%g\t&\t%f", dt_list(i), vy_error(i), order_vy(i));
end
fprintf("\nAverage order: %f", avg_order_vy);

%solve for dt = 0.1 and draw results
for i= 1:length(t_final)
    [~,~,~,~, x, y] = trapezoidal_method(t0, t_final(i), x0, y0, vx0, vy0, dt_list(2), f_x, f_y, f_vx, f_vy, a, b, c, d, r, alpha, beta, g);
    
    %draw results
    figure;
    for i = 1:length(x)
        draw_disk(x(i), y(i), r);
        axis([a b, c d]);
        axis square;
        hold on;
        plot(x,y, '.-');
        hold off;
        pause(dt_list(2));
    end
end

