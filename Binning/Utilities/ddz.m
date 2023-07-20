% function: ddz
% cite: Schneiders, J.F.G. & Scarano, F. Exp Fluids (2016) 57:139. doi:10.1007/s00348-016-2225-6
% http://link.springer.com/article/10.1007/s00348-016-2225-6
 
function dudz = ddz(dz,u,order)
 
if nargin == 2
     
    dudz = zeros(size(u));
     
    % 1st order Upwind on boundaries
     
    dudz(:,:,1)  = 1/1/dz * (u(:,:,2) - u(:,:,1));
    dudz(:,:,end) = 1/1/dz * (u(:,:,end) - u(:,:,end-1));
     
    % 2nd order Central difference on 2nd cells
     
    dudz(:,:,2:end-1) = 1/2/dz * ( - u(:,:,1:end-2) + u(:,:,3:end) );
     
else
     
    rkw = [0        0         0       -1/2  0   1/2     0       0       0
        0        0         1/12    -2/3 0   2/3     -1/12   0       0
        0        -1/60     3/20    -3/4 0   3/4     -3/20   1/60    0
        1/280   -4/105    1/5     -4/5  0   4/5     -1/5    4/105   -1/280];
     
    dudz = zeros(size(u));
     
    [~,~,nz,~] = size(u);
     
    if order == 10
         
        % Upwind
         
        dudz(:,:,1:nz-1) = 1/1/dz * (u(:,:,2:nz) - u(:,:,1:nz-1));
        dudz(:,:,nz)     = 1/1/dz * (u(:,:,nz) - u(:,:,nz-1));
         
    elseif order == 20
         
        % 2nd order Upwind on boundaries
         
        dudz(:,:,1)   = 1/dz * (-3/2*u(:,:,1) + 2*u(:,:,2) - 1/2*u(:,:,3));
        dudz(:,:,nz)  = -1/dz * (-3/2*u(:,:,nz) + 2*u(:,:,nz-1) - 1/2*u(:,:,nz-2));
         
        % 2nd order Central difference on 2nd cells
         
        dudz(:,:,2:nz-1) = 1/dz * ( rkw(1,4) * u(:,:,1:nz-2) + rkw(1,6) * u(:,:,3:nz) );
         
    elseif order == 200
         
        % 1st order Upwind on boundaries
         
        dudz(:,:,1)  = 1/1/dz * (u(:,:,2) - u(:,:,1));
        dudz(:,:,nz) = 1/1/dz * (u(:,:,nz) - u(:,:,nz-1));
         
        % 2nd order Central difference on 2nd cells
         
        dudz(:,:,2:nz-1) = 1/dz * ( rkw(1,4) * u(:,:,1:nz-2) + rkw(1,6) * u(:,:,3:nz) );
         
    elseif order == 40
         
        % 4th order Upwind on boundaries
         
        dudz(:,:,1:2)     = 1/dz * (-25/12*u(:,:,1:2) + 4*u(:,:,2:3) - 3*u(:,:,3:4) + 4/3*u(:,:,4:5) - 1/4*u(:,:,5:6));
        dudz(:,:,nz-1:nz) = -1/dz * (-25/12*u(:,:,nz-1:nz) + 4*u(:,:,nz-2:nz-1) - 3*u(:,:,nz-3:nz-2) + 4/3*u(:,:,nz-4:nz-3) - 1/4*u(:,:,nz-5:nz-4));
         
        % 4nd order Central difference on 3rd cells
         
        dudz(:,:,3:nz-2) = 1/dz * ( rkw(2,3) * u(:,:,1:nz-4) + rkw(2,4) * u(:,:,2:nz-3) ...
            + rkw(2,5) * u(:,:,3:nz-2) + rkw(2,6) * u(:,:,4:nz-1) + rkw(2,7) * u(:,:,5:nz) );
         
    elseif order == 400
         
        % 2nd order Upwind on boundaries
         
        dudz(:,:,1)  = +1/dz * (-3/2*u(:,:,1) + 2*u(:,:,2) - 1/2*u(:,:,3));
        dudz(:,:,nz) = -1/dz * (-3/2*u(:,:,nz) + 2*u(:,:,nz-1) - 1/2*u(:,:,nz-2));
         
        % 2nd order Central difference on 2nd cells
         
        dudz(:,:,[2 nz-1]) = 1/dz * ( rkw(1,4) * u(:,:,[1 nz-2]) + rkw(1,6) * u(:,:,[3 nz]) );
         
        % 4nd order Central difference on 3rd cells
         
        dudz(:,:,3:nz-2) = 1/dz * ( rkw(2,3) * u(:,:,1:nz-4) + rkw(2,4) * u(:,:,2:nz-3) ...
            + rkw(2,5) * u(:,:,3:nz-2) + rkw(2,6) * u(:,:,4:nz-1) + rkw(2,7) * u(:,:,5:nz) );
         
    elseif order == 4000
         
        % 1st order Upwind on boundaries
         
        dudz(:,:,1)  = 1/1/dz * (u(:,:,2) - u(:,:,1));
        dudz(:,:,nz) = 1/1/dz * (u(:,:,nz) - u(:,:,nz-1));
         
        % 2nd order Central difference on 2nd cells
         
        dudz(:,:,[2 nz-1]) = 1/dz * ( rkw(1,4) * u(:,:,[1 nz-2]) + rkw(1,6) * u(:,:,[3 nz]) );
         
        % 4nd order Central difference on 3rd cells
         
        dudz(:,:,3:nz-2) = 1/dz * ( rkw(2,3) * u(:,:,1:nz-4) + rkw(2,4) * u(:,:,2:nz-3) ...
            + rkw(2,5) * u(:,:,3:nz-2) + rkw(2,6) * u(:,:,4:nz-1) + rkw(2,7) * u(:,:,5:nz) );
         
    elseif order == 888
         
        % 1st order Upwind on boundaries
         
        dudz(:,:,1)  = 1/1/dz * (u(:,:,2) - u(:,:,1));
        dudz(:,:,nz) = 1/1/dz * (u(:,:,nz) - u(:,:,nz-1));
         
        % 2nd order Central difference on 2nd cells
         
        dudz(:,:,[2 nz-1]) = 1/dz * ( rkw(1,4) * u(:,:,[1 nz-2]) + rkw(1,6) * u(:,:,[3 nz]) );
         
        % Richard
         
        dudz(:,:,3:nz-2) = 1/12/dz * ( 1 * u(:,:,1:nz-4) - 8 * u(:,:,2:nz-3) ...
            + 0 * u(:,:,3:nz-2) + 8 * u(:,:,4:nz-1) - 1 * u(:,:,5:nz) );
         
    elseif order == 999
         
        % 1st order Upwind on boundaries
         
        dudz(:,:,1)  = 1/1/dz * (u(:,:,2) - u(:,:,1));
        dudz(:,:,nz) = 1/1/dz * (u(:,:,nz) - u(:,:,nz-1));
         
        % 2nd order Central difference on 2nd cells
         
        dudz(:,:,[2 nz-1]) = 1/dz * ( rkw(1,4) * u(:,:,[1 nz-2]) + rkw(1,6) * u(:,:,[3 nz]) );
         
        % Least Squares
         
        dudz(:,:,3:nz-2) = 1/10/dz * ( -2 * u(:,:,1:nz-4) - 1 * u(:,:,2:nz-3) ...
            + 0 * u(:,:,3:nz-2) + 1 * u(:,:,4:nz-1) + 2 * u(:,:,5:nz) );
                 
    elseif order == 60
         
        % Upwind on boundaries
         
        dudz(:,:,1)       = 1/1/dz * (u(:,:,2) - u(:,:,1));
        dudz(:,:,nz)     = 1/1/dz * (u(:,:,nz) - u(:,:,nz-1));
         
        % 2nd order Central difference on 2nd cells
         
        dudz(:,:,[2 nz-1]) = 1/dz * ( rkw(1,4) * u(:,:,[1 nz-2]) + rkw(1,6) * u(:,:,[3 nz]) );
         
        % 4nd order Central difference on 3rd cells
         
        dudz(:,:,[3 nz-2]) = 1/dz * ( rkw(2,3) * u(:,:,[1 nz-4]) + rkw(2,4) * u(:,:,[2 nz-3]) ...
            + rkw(2,5) * u(:,:,[3 nz-2]) + rkw(2,6) * u(:,:,[4 nz-1]) + rkw(2,7) * u(:,:,[5 nz]) );
         
        % 6th order Central difference on interior cells
         
        dudz(:,:,4:nz-3) = 1/dz * ( rkw(3,2) * u(:,:,1:nz-6) + rkw(3,3) * u(:,:,2:nz-5) + rkw(3,4) * u(:,:,3:nz-4) ...
            + rkw(3,5) * u(:,:,4:nz-3) + rkw(3,6) * u(:,:,5:nz-2) + rkw(3,7) * u(:,:,6:nz-1) + rkw(3,8) * u(:,:,7:nz) );
         
    elseif order == 80
         
        % Upwind on boundaries
         
        dudz(:,:,1)       = 1/1/dz * (u(:,:,2) - u(:,:,1));
        dudz(:,:,nz)     = 1/1/dz * (u(:,:,nz) - u(:,:,nz-1));
         
        % 2nd order Central difference on 2nd cells
         
        dudz(:,:,[2 nz-1]) = 1/dz * ( rkw(1,4) * u(:,:,[1 nz-2]) + rkw(1,6) * u(:,:,[3 nz]) );
         
        % 4nd order Central difference on 3rd cells
         
        dudz(:,:,[3 nz-2]) = 1/dz * ( rkw(2,3) * u(:,:,[1 nz-4]) + rkw(2,4) * u(:,:,[2 nz-3]) ...
            + rkw(2,5) * u(:,:,[3 nz-2]) + rkw(2,6) * u(:,:,[4 nz-1]) + rkw(2,7) * u(:,:,[5 nz]) );
         
        % 6th order Central difference on interior cells
         
        dudz(:,:,[4 nz-3]) = 1/dz * ( rkw(3,2) * u(:,:,[1 nz-6]) + rkw(3,3) * u(:,:,[2 nz-5]) + rkw(3,4) * u(:,:,[3 nz-4]) ...
            + rkw(3,5) * u(:,:,[4 nz-3]) + rkw(3,6) * u(:,:,[5 nz-2]) + rkw(3,7) * u(:,:,[6 nz-1]) + rkw(3,8) * u(:,:,[7 nz]) );
         
        % 8th order Central difference on interior cells
         
        dudz(:,:,5:nz-4) = 1/dz * ( ...
            rkw(4,1) * u(:,:,1:nz-8) + ...
            rkw(4,2) * u(:,:,2:nz-7) + ...
            rkw(4,3) * u(:,:,3:nz-6) + ...
            rkw(4,4) * u(:,:,4:nz-5) + ...
            rkw(4,5) * u(:,:,5:nz-4) + ...
            rkw(4,6) * u(:,:,6:nz-3) + ...
            rkw(4,7) * u(:,:,7:nz-2) + ...
            rkw(4,8) * u(:,:,8:nz-1) + ...
            rkw(4,9) * u(:,:,9:nz) );
         
    elseif order == 15
         
        % 1st order Upwind on boundaries
         
        dudz(:,:,1)  = 1/1/dz * (u(:,:,2) - u(:,:,1));
        dudz(:,:,nz) = 1/1/dz * (u(:,:,nz) - u(:,:,nz-1));
         
        % 2nd order Central difference on 2nd cells
         
        dudz(:,:,[2 nz-1]) = 1/dz * ( rkw(1,4) * u(:,:,[1 nz-2]) + rkw(1,6) * u(:,:,[3 nz]) );
         
        % 4nd order Central difference on 3rd cells
         
        dudz(:,:,3:nz-2) = 1/dz/10 * ( -2 * u(:,:,1:nz-4)  - 1 * u(:,:,2:nz-3) ...
            + 1 * u(:,:,4:nz-1) + 2 * u(:,:,5:nz) );
         
    elseif order == 123
         
        % 1st order Upwind on boundaries
         
        dudz(:,:,1)  = 1/1/dz * (u(:,:,2) - u(:,:,1));
        dudz(:,:,nz) = 1/1/dz * (u(:,:,nz) - u(:,:,nz-1));
         
        % 2nd order Central difference on 2nd cells
         
        dudz(:,:,2:nz-1) = 1/dz * ( rkw(1,4) * u(:,:,1:nz-2) + rkw(1,6) * u(:,:,3:nz) );
         
        % CR4s on interior cells
         
        A  = 1239;
        A1 = 272;
        A2 = 1036;
        A4 = 0;
        A8 = -69;
         
        dudz(:,:,9:nz-8) = 1/A * ( ...
            A1/2/1/dz * (u(:,:,9+1:nz-8+1) - u(:,:,9-1:nz-8-1)) + ...
            A2/2/2/dz * (u(:,:,9+2:nz-8+2) - u(:,:,9-2:nz-8-2)) + ...
            A4/2/4/dz * (u(:,:,9+4:nz-8+4) - u(:,:,9-4:nz-8-4)) + ...
            A8/2/8/dz * (u(:,:,9+8:nz-8+8) - u(:,:,9-8:nz-8-8)) ...
            );
         
    else
         
        error('select derivative scheme')
         
         
    end
end