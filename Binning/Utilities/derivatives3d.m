function [dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz] ...
    = derivatives3d(hx, hy, hz, u, dm)
 
if nargin == 4, dm = 200; end
 
dudx = ddx(hx, u(:,:,:,1), dm);
dvdx = ddx(hx, u(:,:,:,2), dm);
dwdx = ddx(hx, u(:,:,:,3), dm);
 
dudy = ddy(hy, u(:,:,:,1), dm);
dvdy = ddy(hy, u(:,:,:,2), dm);
dwdy = ddy(hy, u(:,:,:,3), dm);
 
dudz = ddz(hz, u(:,:,:,1), dm);
dvdz = ddz(hz, u(:,:,:,2), dm);
dwdz = ddz(hz, u(:,:,:,3), dm);