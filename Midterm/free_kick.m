function out = free_kick(t, y)
%parameters
g = -9.81;
m = 0.437;
c_drag = 0.0057;
c_lift = 0.0061;

%assign temporary variables to preserve my sanity
p_x = y(1);
p_y = y(2);
p_z = y(3);
v_x = y(4);
v_y = y(5);
v_z = y(6);
s_x = y(7);
s_y = y(8);
s_z = y(9);

%calculate the magnitude of the velocity vector
v_mag = sqrt(v_x^2 + v_y^2 + v_z^2);

out(1) = v_x;
out(2) = v_y;
out(3) = v_z;
out(4) = -c_drag*v_mag*v_x/m + c_lift*v_mag*(s_y*v_z - s_z*v_y)/m;
out(5) = -c_drag*v_mag*v_y/m + c_lift*v_mag*(s_z*v_x - s_x*v_z)/m;
out(6) = g - c_drag*v_mag*v_z/m + c_lift*v_mag*(s_x*v_y - s_y*v_x)/m;
out(7) = 0;
out(8) = 0;
out(9) = 0;
end



