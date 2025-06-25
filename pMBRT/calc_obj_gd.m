function [y1] = calc_obj_gd(d,wp)
w_i = wp.r_i;
w_x = wp.r_x;
w_z = wp.r_z;

intp = wp.intp;
tvx=wp.tvx;
tvz=wp.tvz;

Ixd = intp*d;
y1 = 0;
y1 = y1+(w_i/size(intp,1))*norm(Ixd,1);
y1 = y1+(w_x/size(tvx,1))*norm(tvx*Ixd,1);
y1 = y1+(w_z/size(tvz,1))*norm(tvz*Ixd,1);
end