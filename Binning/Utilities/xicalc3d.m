%% Function to calculate vorticity field in 3d
 
function xi = xicalc3d(hx,hy,hz,u,dm)
 
if nargin == 4, dm = 200; end
 
[~, dudy, dudz, dvdx, ~, dvdz, dwdx, dwdy, ~] = derivatives3d(hx, hy, hz, u, dm);
 
xi(:,:,:,1) = dwdy - dvdz;
xi(:,:,:,2) = dudz - dwdx;
xi(:,:,:,3) = dvdx - dudy;