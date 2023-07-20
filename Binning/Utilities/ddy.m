% function: ddy
% cite: Schneiders, J.F.G. & Scarano, F. Exp Fluids (2016) 57:139. doi:10.1007/s00348-016-2225-6
% http://link.springer.com/article/10.1007/s00348-016-2225-6

function dudy = ddy(dy,u,order)

if nargin == 2
    dudy = zeros(size(u));
    
    % 1st order Upwind on boundaries
    
    dudy(1,:,:)  = 1/1/dy * (u(2,:,:) - u(1,:,:));
    dudy(end,:,:) = 1/1/dy * (u(end,:,:) - u(end-1,:,:));
    
    % 2nd order Central difference on 2nd cells
    
    dudy(2:end-1,:,:) = 1/2/dy * ( - u(1:end-2,:,:) + u(3:end,:,:) );
    
    
else
    
    rkw = [0        0         0       -1/2  0   1/2     0       0       0
        0        0         1/12    -2/3 0   2/3     -1/12   0       0
        0        -1/60     3/20    -3/4 0   3/4     -3/20   1/60    0
        1/280   -4/105    1/5     -4/5  0   4/5     -1/5    4/105   -1/280];
    
    dudy = zeros(size(u));
    
    [ny,~,~,~] = size(u);
    
    if order == 10
        
        % Upwind
        
        dudy(1:ny-1,:,:) = 1/1/dy * (u(2:ny,:,:) - u(1:ny-1,:,:));
        dudy(ny,:,:)     = 1/1/dy * (u(ny,:,:) - u(ny-1,:,:));
        
    elseif order == 20
        
        % Upwind on boundaries
        
        dudy(1,:,:)   = 1/dy * (-3/2*u(1,:,:) + 2*u(2,:,:) - 1/2*u(3,:,:));
        dudy(ny,:,:)  = -1/dy * (-3/2*u(ny,:,:) + 2*u(ny-1,:,:) - 1/2*u(ny-2,:,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudy(2:ny-1,:,:) = 1/dy * ( rkw(1,4) * u(1:ny-2,:,:) + rkw(1,6) * u(3:ny,:,:) );
        
    elseif order == 200
        
        % 1st order Upwind on boundaries
        
        dudy(1,:,:)  = 1/1/dy * (u(2,:,:) - u(1,:,:));
        dudy(ny,:,:) = 1/1/dy * (u(ny,:,:) - u(ny-1,:,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudy(2:ny-1,:,:) = 1/dy * ( rkw(1,4) * u(1:ny-2,:,:) + rkw(1,6) * u(3:ny,:,:) );
        
    elseif order == 40
        
        % 4th order Upwind on boundaries
        
        dudy(1:2,:,:)     = 1/dy * (-25/12*u(1:2,:,:) + 4*u(2:3,:,:) - 3*u(3:4,:,:) + 4/3*u(4:5,:,:) - 1/4*u(5:6,:,:));
        dudy(ny-1:ny,:,:) = -1/dy * (-25/12*u(ny-1:ny,:,:) + 4*u(ny-2:ny-1,:,:) - 3*u(ny-3:ny-2,:,:) + 4/3*u(ny-4:ny-3,:,:) - 1/4*u(ny-5:ny-4,:,:));
        
        % 4nd order Central difference on 3rd cells
        
        dudy(3:ny-2,:,:) = 1/dy * ( rkw(2,3) * u(1:ny-4,:,:) + rkw(2,4) * u(2:ny-3,:,:) ...
            + rkw(2,5) * u(3:ny-2,:,:) + rkw(2,6) * u(4:ny-1,:,:) + rkw(2,7) * u(5:ny,:,:) );
        
    elseif order == 400
        
        % 2nd order Upwind on boundaries
        
        dudy(1,:,:)  = +1/dy * (-3/2*u(1,:,:) + 2*u(2,:,:) - 1/2*u(3,:,:));
        dudy(ny,:,:) = -1/dy * (-3/2*u(ny,:,:) + 2*u(ny-1,:,:) - 1/2*u(ny-2,:,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudy([2 ny-1],:,:) = 1/dy * ( rkw(1,4) * u([1 ny-2],:,:) + rkw(1,6) * u([3 ny],:,:) );
        
        % 4nd order Central difference on 3rd cells
        
        dudy(3:ny-2,:,:) = 1/dy * ( rkw(2,3) * u(1:ny-4,:,:) + rkw(2,4) * u(2:ny-3,:,:) ...
            + rkw(2,5) * u(3:ny-2,:,:) + rkw(2,6) * u(4:ny-1,:,:) + rkw(2,7) * u(5:ny,:,:) );
        
    elseif order == 4000
        
        % 1st order Upwind on boundaries
        
        dudy(1,:,:)  = 1/1/dy * (u(2,:,:) - u(1,:,:));
        dudy(ny,:,:) = 1/1/dy * (u(ny,:,:) - u(ny-1,:,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudy([2 ny-1],:,:) = 1/dy * ( rkw(1,4) * u([1 ny-2],:,:) + rkw(1,6) * u([3 ny],:,:) );
        
        % 4nd order Central difference on 3rd cells
        
        dudy(3:ny-2,:,:) = 1/dy * ( rkw(2,3) * u(1:ny-4,:,:) + rkw(2,4) * u(2:ny-3,:,:) ...
            + rkw(2,5) * u(3:ny-2,:,:) + rkw(2,6) * u(4:ny-1,:,:) + rkw(2,7) * u(5:ny,:,:) );
        
    elseif order == 888
        
        % 1st order Upwind on boundaries
        
        dudy(1,:,:)  = 1/1/dy * (u(2,:,:) - u(1,:,:));
        dudy(ny,:,:) = 1/1/dy * (u(ny,:,:) - u(ny-1,:,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudy([2 ny-1],:,:) = 1/dy * ( rkw(1,4) * u([1 ny-2],:,:) + rkw(1,6) * u([3 ny],:,:) );
        
        % Richardson
        
        dudy(3:ny-2,:,:) = 1/12/dy * ( 1 * u(1:ny-4,:,:) - 8 * u(2:ny-3,:,:) ...
            + 0 * u(3:ny-2,:,:) + 8 * u(4:ny-1,:,:) - 1 * u(5:ny,:,:) );
        
    elseif order == 999
        
        % 1st order Upwind on boundaries
        
        dudy(1,:,:)  = 1/1/dy * (u(2,:,:) - u(1,:,:));
        dudy(ny,:,:) = 1/1/dy * (u(ny,:,:) - u(ny-1,:,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudy([2 ny-1],:,:) = 1/dy * ( rkw(1,4) * u([1 ny-2],:,:) + rkw(1,6) * u([3 ny],:,:) );
        
        % Least Squares
        
        dudy(3:ny-2,:,:) = 1/10/dy * ( -2 * u(1:ny-4,:,:) -1 * u(2:ny-3,:,:) ...
            + 0 * u(3:ny-2,:,:) + 1 * u(4:ny-1,:,:) + 2 * u(5:ny,:,:) );
        
    elseif order == 4000
        
        % 1st order Upwind on boundaries
        
        dudy(1,:,:)  = 1/1/dy * (u(2,:,:) - u(1,:,:));
        dudy(ny,:,:) = 1/1/dy * (u(ny,:,:) - u(ny-1,:,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudy([2 ny-1],:,:) = 1/dy * ( rkw(1,4) * u([1 ny-2],:,:) + rkw(1,6) * u([3 ny],:,:) );
        
        % 4nd order Central difference on 3rd cells
        
        dudy(3:ny-2,:,:) = 1/dy * ( rkw(2,3) * u(1:ny-4,:,:) + rkw(2,4) * u(2:ny-3,:,:) ...
            + rkw(2,5) * u(3:ny-2,:,:) + rkw(2,6) * u(4:ny-1,:,:) + rkw(2,7) * u(5:ny,:,:) );
        
    elseif order == 60
        
        % Upwind on boundaries
        
        dudy(1,:,:)       = 1/1/dy * (u(2,:,:) - u(1,:,:));
        dudy(ny,:,:)     = 1/1/dy * (u(ny,:,:) - u(ny-1,:,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudy([2 ny-1],:,:) = 1/dy * ( rkw(1,4) * u([1 ny-2],:,:) + rkw(1,6) * u([3 ny],:,:) );
        
        % 4nd order Central difference on 3rd cells
        
        dudy([3 ny-2],:,:) = 1/dy * ( rkw(2,3) * u([1 ny-4],:,:) + rkw(2,4) * u([2 ny-3],:,:) ...
            + rkw(2,5) * u([3 ny-2],:,:) + rkw(2,6) * u([4 ny-1],:,:) + rkw(2,7) * u([5 ny],:,:) );
        
        % 6th order Central difference on interior cells
        
        dudy(4:ny-3,:,:) = 1/dy * ( rkw(3,2) * u(1:ny-6,:,:) + rkw(3,3) * u(2:ny-5,:,:) + rkw(3,4) * u(3:ny-4,:,:) ...
            + rkw(3,5) * u(4:ny-3,:,:) + rkw(3,6) * u(5:ny-2,:,:) + rkw(3,7) * u(6:ny-1,:,:) + rkw(3,8) * u(7:ny,:,:) );
        
    elseif order == 80
        
        % Upwind on boundaries
        
        dudy(1,:,:)       = 1/1/dy * (u(2,:,:) - u(1,:,:));
        dudy(ny,:,:)     = 1/1/dy * (u(ny,:,:) - u(ny-1,:,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudy([2 ny-1],:,:) = 1/dy * ( rkw(1,4) * u([1 ny-2],:,:) + rkw(1,6) * u([3 ny],:,:) );
        
        % 4nd order Central difference on 3rd cells
        
        dudy([3 ny-2],:,:) = 1/dy * ( rkw(2,3) * u([1 ny-4],:,:) + rkw(2,4) * u([2 ny-3],:,:) ...
            + rkw(2,5) * u([3 ny-2],:,:) + rkw(2,6) * u([4 ny-1],:,:) + rkw(2,7) * u([5 ny],:,:) );
        
        % 6th order Central difference on interior cells
        
        dudy([4 ny-3],:,:) = 1/dy * ( rkw(3,2) * u([1 ny-6],:,:) + rkw(3,3) * u([2 ny-5],:,:) + rkw(3,4) * u([3 ny-4],:,:) ...
            + rkw(3,5) * u([4 ny-3],:,:) + rkw(3,6) * u([5 ny-2],:,:) + rkw(3,7) * u([6 ny-1],:,:) + rkw(3,8) * u([7 ny],:,:) );
        
        % 8th order Central difference on interior cells
        
        dudy(5:ny-4,:,:) = 1/dy * ( ...
            rkw(4,1) * u(1:ny-8,:,:) + ...
            rkw(4,2) * u(2:ny-7,:,:) + ...
            rkw(4,3) * u(3:ny-6,:,:) + ...
            rkw(4,4) * u(4:ny-5,:,:) + ...
            rkw(4,5) * u(5:ny-4,:,:) + ...
            rkw(4,6) * u(6:ny-3,:,:) + ...
            rkw(4,7) * u(7:ny-2,:,:) + ...
            rkw(4,8) * u(8:ny-1,:,:) + ...
            rkw(4,9) * u(9:ny,:,:) );
        
    elseif order == 15
        
        % 1st order Upwind on boundaries
        
        dudy(1,:,:)  = 1/1/dy * (u(2,:,:) - u(1,:,:));
        dudy(ny,:,:) = 1/1/dy * (u(ny,:,:) - u(ny-1,:,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudy([2 ny-1],:,:) = 1/dy * ( rkw(1,4) * u([1 ny-2],:,:) + rkw(1,6) * u([3 ny],:,:) );
        
        % 4nd order Central difference on 3rd cells
        
        dudy(3:ny-2,:,:) = 1/dy/10 * ( -2 * u(1:ny-4,:,:) -1 * u(2:ny-3,:,:) ...
            + 1 * u(4:ny-1,:,:) + 2 * u(5:ny,:,:) );
        
    elseif order == 123
        
        % 1st order Upwind on boundaries
        
        dudy(1,:,:)  = 1/1/dy * (u(2,:,:) - u(1,:,:));
        dudy(ny,:,:) = 1/1/dy * (u(ny,:,:) - u(ny-1,:,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudy(2:ny-1,:,:) = 1/dy * ( rkw(1,4) * u(1:ny-2,:,:) + rkw(1,6) * u(3:ny,:,:) );
        
        % CR4s on interior cells
        
        A  = 1239;
        A1 = 272;
        A2 = 1036;
        A4 = 0;
        A8 = -69;
        
        dudy(9:ny-8,:,:) = 1/A * ( ...
            A1/2/1/dy * (u(9+1:ny-8+1,:,:) - u(9-1:ny-8-1,:,:)) + ...
            A2/2/2/dy * (u(9+2:ny-8+2,:,:) - u(9-2:ny-8-2,:,:)) + ...
            A4/2/4/dy * (u(9+4:ny-8+4,:,:) - u(9-4:ny-8-4,:,:)) + ...
            A8/2/8/dy * (u(9+8:ny-8+8,:,:) - u(9-8:ny-8-8,:,:)) ...
            );
        
    else
        
        error('select derivative scheme')
        
    end
    
    
    
    
end