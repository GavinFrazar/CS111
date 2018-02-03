close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the main file for part 4 of the midterm assignment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters
H_wall = 2;
L_goal = 7.32;
H_goal = 2.44;

%initial conditions
t_start = 0;
t_final = 2;
dt = 0.02;
x0 = [-15 25.5 -4 16];
y0 = [23 15 35 28];
z0 = [0 0 0 0];
vx0 = [18 -28 -8 -25];
vy0 = [-23 -11 -36 -20];
vz0 = [8.5 7.5 5 8];
sx0 = [0.10 0 -0.32 0.10];
sy0 = [0.10 0 0 0.15];
sz0 = [-0.99 1 0.95 0.98];
xw_start = [-11.2 18.1 -4.3 12.6];
yw_start = [14.7 9.7 25.9 19.6];
xw_end = [-8.6 17 -0.8 9.9];
yw_end = [16 11.3 25.9 20.8];

%combine initial conditions into matrices
IC = vertcat(x0,y0,z0,vx0,vy0,vz0,sx0,sy0,sz0);
wall = vertcat(xw_start, yw_start, xw_end, yw_end);

%returns number of columns in a matrix
cols =@(x) size(x,2);

%find trajectory of the ball using RK4 (before checking for collisions)
for c= 1:cols(IC)
    [t,sol] = solve_system_RK4(@free_kick, IC(:,c), t_start, t_final, dt);
    
    %analyze the trajectory to check for collisions
    [hit_wall, hit_goal] = analyze_trajectory(t, sol, wall(:,c), dt);
    
    %plot results
    figure;
    plot_trajectory(sol(1,:),sol(2,:),sol(3,:), xw_start(c), yw_start(c), xw_end(c), yw_end(c), H_wall);
    
    %print results of analysis
    fprintf("\nTrial: %u",c);
    if (hit_wall)
       fprintf("\nBall hit the wall.\n"); 
    elseif (hit_goal)
            fprintf("\nGoal!\n");
    else
            fprintf("\nMissed both wall and goal.\n");
    end
end