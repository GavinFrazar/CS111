function [c_new] = upwind_scheme(t_start,t_final,C,x,y,Nx,Ny,f,boundary_condition,vel,c_start)
    dx = x(2)-x(1);
    dy = y(2)-y(1);  
      
    %impose initial conditions
    c_old = c_start(x,y);
    c_new = c_old; %for the case where t_start == t_final and nothing changes
    t = t_start; 
    [vx,vy] = vel(t,x,y); %in our cases velocity doesn't actually depend on
                          %time which we can exploit to speed things up by 
                          %calculating velocity at each (x,y) once.
                          %This will fail if velocity depends on time, in
                          %which case we just have to move this statement
                          %and the logical masking into the time-loop to
                          %recalculate at each new time t.

    dt = C*dx*dy./(vy*dx+vx*dy); %dt is now a matrix of all possible dt
    idx = dt > 0;                %this mask flags all dt > 0 in the matrix
    dt = min(min(abs(idx.*dt + ~idx.*2)));   %we now get the maximum dt which still 
                                        %satisfies the stability condition everywhere.
    
    %initialize masking matrices which will be borders of an Nx*Ny matrix
    LEFT = zeros(Nx,Ny);
    RIGHT = LEFT;
    BOTTOM = LEFT;
    TOP = LEFT;
    
    %set column vectors for the left and right of an Nx*Ny matrix
    LEFT(:,1) = 1;
    RIGHT(:,end) = 1;
    
    %set row vectors for the bottom and top of an Nx*Ny matrix
    BOTTOM(end,:) = 1;
    TOP(1,:) = 1;
    
    %logical matrix masks
    idx1 = (vx >= 0) & ~TOP; %use backward difference in x-axis at these points. 
    idx2 = (vx < 0) & ~BOTTOM; %use forward difference in x-axis.
    idx3 = ((vy >=0) & ~LEFT); %use backward difference in y-axis.
    idx4 = ((vy < 0) & ~RIGHT); %use forward difference in y-axis.
    idx5 = ~(idx1 | idx2);  %use boundary conditions where neither
                            %the forward nor backward difference will work
                            %(x-axis)
    idx6 = ~(idx3 | idx4);  %use boundary conditions (y-axis)
    
    %use zero-padding and Matlab's diff() function to apply the forward and 
    %backward difference. The results will be inaccurate at the boundary of
    %course, but the logical masks will prevent those boundary results from
    %being applied in the final calculation.
    row_pad = zeros(1,size(c_old,2));
    col_pad = zeros(size(c_old,1),1);
    while t < t_final
        %ensure we step into t_final exactly
        if (t+dt > t_final)
            dt = t_final - t;
        end
        %matrix operation to find derivative of c(t,x,y) w/ respect to x
        %at all grid points.
        dcdx = (idx1.*(diff([row_pad; c_old],1,1)))/dx...
               + (idx2.*(diff([c_old; row_pad],1,1)))/dx...
               + idx5.*(boundary_condition(t,x,y));
            
        %matrix operation to find derivative of c(t,x,y) w/ respect to y at
        %all grid points.
        dcdy = (idx3.*(diff([col_pad, c_old],1,2)))/dy...
            + (idx4.*(diff([c_old, col_pad],1,2)))/dy...
            + idx6.*(boundary_condition(t,x,y));
            
        c_new = c_old - dt*(vx.*dcdx + vy.*dcdy - f(t,x,y));
        c_old = c_new;
        t = t + dt;
    end
end