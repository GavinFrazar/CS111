function plot_trajectory(x,y,z, wall_x_start, wall_y_start, wall_x_end, wall_y_end, wall_height)

n_cut = length(x);
n_cut = min([n_cut, min(find(y<0))]);
n_cut = min([n_cut, min(find(z<-1.e-12))]);

x_field = [ -40, -40, 40, 40 ];
y_field = [ 40, 0, 0, 40 ];
z_field = [ 0, 0, 0, 0];

x_pen_box = [-20, -20, 20, 20];
y_pen_box = [0, 16.5, 16.5, 0];
z_pen_box = [0, 0, 0, 0];

x_goal_box = [-9.16, -9.16, 9.16, 9.16];
y_goal_box = [0, 5.5, 5.5, 0];
z_goal_box = [0, 0, 0, 0];

x_goal = [-3.66, -3.66, 3.66, 3.66];
y_goal = [ 0, 0, 0, 0];
z_goal = [0, 2.44, 2.44, 0];

x_wall = [wall_x_start, wall_x_start, wall_x_end, wall_x_end];
y_wall = [wall_y_start, wall_y_start, wall_y_end, wall_y_end];
z_wall = [0, wall_height, wall_height, 0];

alpha = asin(5.5/9.15);
theta = linspace(alpha, pi-alpha, 20);
x_arc =  0*ones(size(theta)) + 9.15*cos(theta);
y_arc = 11*ones(size(theta)) + 9.15*sin(theta);
z_arc = zeros(size(theta));

theta = linspace(0, 2*pi, 20);
x_pen =  0*ones(size(theta)) + 0.15*cos(theta);
y_pen = 11*ones(size(theta)) + 0.15*sin(theta);
z_pen = zeros(size(theta));

subplot(2,2,[1 2]);
plot3(x(1:n_cut), y(1:n_cut), z(1:n_cut), '.-');
hold on
fill3(x_field, y_field, z_field, 'white');
plot3(x_field, y_field, z_field, 'k', 'LineWidth', 2);
plot3(x_pen_box, y_pen_box, z_pen_box, 'k', 'LineWidth', 2);
plot3(x_goal_box, y_goal_box, z_goal_box, 'k', 'LineWidth', 2);
plot3(x_arc, y_arc, z_arc, 'k', 'LineWidth', 2);
plot3(x_pen, y_pen, z_pen, 'k', 'LineWidth', 2);
plot3(x_goal, y_goal, z_goal, 'k', 'LineWidth', 2);
fill3(x_goal, y_goal, z_goal, 'green');
plot3(x_wall, y_wall, z_wall, 'k', 'LineWidth', 2);
fill3(x_wall, y_wall, z_wall, 'red');
hold off

axis ([-40, 40, 0, 40, 0, 4]);
axis equal;
% axis off;

% ax.View = [0, 90];
set(gca,'View',[180-atan2(x(1),y(1))*180/pi, 20]);
camproj('perspective')
% camzoom(2)

subplot(2,2,3);
plot(x_field, y_field, 'k', 'LineWidth', 2);
title('Top view');
hold on
plot(x_pen_box, y_pen_box, 'k', 'LineWidth', 2);
plot(x_goal_box, y_goal_box, 'k', 'LineWidth', 2);
plot(x_arc, y_arc, 'k', 'LineWidth', 2);
plot(x_pen, y_pen, 'k', 'LineWidth', 2);
plot(x_goal, y_goal, 'g', 'LineWidth', 2);
plot(x_wall, y_wall, 'r', 'LineWidth', 2);
plot(x(1:n_cut), y(1:n_cut), '.-');
plot([-3.66, -3.66], [-0.5, 0.5], '-k', 'LineWidth', 2);
plot([3.66, 3.66], [-0.5, 0.5], '-k', 'LineWidth', 2);
hold off

axis ([-40, 40, 0, 40]);
axis equal;
% axis off;


subplot(2,2,4);
plot(y_field, z_field, 'k', 'LineWidth', 2);
title('Side view');
hold on
plot(y_pen_box, z_pen_box, 'k', 'LineWidth', 2);
plot(y_goal_box, z_goal_box, 'k', 'LineWidth', 2);
plot(y_arc, z_arc, 'k', 'LineWidth', 2);
plot(y_goal, z_goal, 'g', 'LineWidth', 3);
plot(y_wall, z_wall, 'k', 'LineWidth', 2);
fill(y_wall, z_wall, 'red');
plot(y(1:n_cut), z(1:n_cut), '.-');
hold off

axis ([0, 40, 0, 4]);
axis equal;
% axis off;

end