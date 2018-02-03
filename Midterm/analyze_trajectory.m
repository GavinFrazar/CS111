function [hit_wall, hit_goal] = analyze_trajectory(t, solutions, wall, dt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function analyzes the trajectory of the ball (part 5 of the  %
%assignment) and determines if there was a collision with the wall %
%or if a goal was scored.                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters
H_wall = 2;
L_goal = 7.32;
H_goal = 2.44;
hit_goal = false;
hit_wall = false;

%macros
x =@(i) solutions(1,i);
y =@(i) solutions(2,i);
z =@(i) solutions(3,i);
xw_start = wall(1);
yw_start = wall(2);
xw_end = wall(3);
yw_end = wall(4);

%calculate slope and y-intercept of the hyperplane that the wall lies on.
m = (yw_end - yw_start)/(xw_end - xw_start);
b = -m*xw_start + yw_start;

for n= 1:(length(t)-1)
    %check for wall collision
    if ( (y(n) > m*x(n)+b) && (y(n+1) < m*x(n+1)+b) )
        %get equations of the hyperplanes that the wall and ball trajectory
        %lie on
        [A_1, B_1, C_1] = get_eq(xw_start, yw_start, xw_end, yw_end);
        [A_2, B_2, C_2] = get_eq(x(n), y(n), x(n+1), y(n+1));
        intersect = find_intersection(A_1,B_1,C_1,A_2,B_2,C_2);
        x_int = intersect(1);
        
        %if x is in the interval of the wall, then there is a potential
        %collision.
        if (x_int >= min(xw_start,xw_end) && x_int <= max(xw_start,xw_end))
            %calculate a smaller dt to step into the point of xy
            %intersection
            dt_temp = dt*(abs(x_int-x(n)))/abs(x(n+1)-x(n));
            
            %get an approximate solution at the moment of xy intersection
            [~, y_temp] = solve_system_RK4(@free_kick, solutions(:,n), t(n), t(n)+dt_temp, dt_temp);
            
            %If z is less than the height of the wall at the point of
            %intersection, then the ball must have collided with the wall.
            if (y_temp(3,2) <= H_wall)
               hit_wall = true;
               hit_goal = false;
               return;
            end
        end        
    end
    
    %check for goal
    if ( y(n) >= 0 && y(n+1) <= 0)
        %get equations of the baseline and the ball trajectory
        [A_1, B_1, C_1] = get_eq(-L_goal/2, 0, L_goal/2, 0);
        [A_2, B_2, C_2] = get_eq(x(n),y(n),x(n+1),y(n+1));
        intersect = find_intersection(A_1,B_1,C_1,A_2,B_2,C_2);
        x_int = intersect(1);
        
        if (abs(x_int) <= L_goal/2)
            dt_temp = dt*(abs(x_int-x(n)))/abs(x(n+1)-x(n));
            [~, y_temp] = solve_system_RK4(@free_kick, solutions(:,n), t(n), t(n)+dt_temp, dt_temp);

            %If z is less than the height of the goal at the point of
            %intersection, then the ball must have gone in the goal.
            if (y_temp(3,2) <= H_goal)
               hit_goal = true;
               return;
            end
        else
            return;
        end
    end
end

end

%get coefficients of the equation of the hyperplane along which the
%two points lie
function [A, B, C] = get_eq(x1, y1, x2, y2)
    A = y2 - y1;
    B = x1 - x2;
    C = A*x1 + B*y1;
end

function intersect = find_intersection(A_1,B_1,C_1,A_2,B_2,C_2)
        %create the matrices of the matrix equation Mx = C
        A = vertcat(A_1, A_2);
        B = vertcat(B_1, B_2);
        C = vertcat(C_1, C_2);
        M = horzcat(A,B);
        
        %solve the matrix equation to obtain the intersection coordinates
        intersect = M\C;
end
