% function: ddx
% cite: Schneiders, J.F.G. & Scarano, F. Exp Fluids (2016) 57:139. doi:10.1007/s00348-016-2225-6
% http://link.springer.com/article/10.1007/s00348-016-2225-6

function dudx = ddx(dx, u, order)

if nargin == 2
    
    dudx = zeros(size(u));
    
    % 1st order Upwind on boundaries
    
    dudx(:,1,:)  = 1/1/dx * (u(:,2,:) - u(:,1,:));
    dudx(:,end,:) = 1/1/dx * (u(:,end,:) - u(:,end-1,:));
    
    % 2nd order Central difference on 2nd cells
    
    dudx(:,2:end-1,:) = 1/2/dx * ( - u(:,1:end-2,:) + u(:,3:end,:) );
    
else
    
    
    rkw = [0        0         0       -1/2  0   1/2     0       0       0
        0        0         1/12    -2/3 0   2/3     -1/12   0       0
        0        -1/60     3/20    -3/4 0   3/4     -3/20   1/60    0
        1/280   -4/105    1/5     -4/5  0   4/5     -1/5    4/105   -1/280];
    
    upw = [-1   1   0       0       0       0   0;
        -3/2    2   -1/2    0       0       0   0;
        -11/6   3   -3/2    1/3     0       0   0;
        -25/12  4   -3      4/3     -1/4    0   0;
        -137/60 5   -5      10/3    -5/4    1/5 0;
        -49/20  6   -15/2   20/3    -15/4   6/5 -1/6];
    
    dudx = zeros(size(u));
    
    [~,nx,~,~] = size(u);
    
    if order == 10
        
        % Upwind
        
        dudx(:,1:nx-1,:) = 1/1/dx * (u(:,2:nx,:) - u(:,1:nx-1,:));
        dudx(:,nx,:)     = 1/1/dx * (u(:,nx,:) - u(:,nx-1,:));
        
    elseif order == 20
        
        % 2nd order Upwind on boundaries
        
        dudx(:,1,:)   = 1/dx * (-3/2*u(:,1,:) + 2*u(:,2,:) - 1/2*u(:,3,:));
        dudx(:,nx,:)  = -1/dx * (-3/2*u(:,nx,:) + 2*u(:,nx-1,:) - 1/2*u(:,nx-2,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudx(:,2:nx-1,:) = 1/dx * ( rkw(1,4) * u(:,1:nx-2,:) + rkw(1,6) * u(:,3:nx,:) );
        
    elseif order == 200
        
        % 1st order Upwind on boundaries
        
        dudx(:,1,:)  = 1/1/dx * (u(:,2,:) - u(:,1,:));
        dudx(:,nx,:) = 1/1/dx * (u(:,nx,:) - u(:,nx-1,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudx(:,2:nx-1,:) = 1/dx * ( rkw(1,4) * u(:,1:nx-2,:) + rkw(1,6) * u(:,3:nx,:) );
        
    elseif order == 40
        
        % 4th order Upwind on boundaries
        
        dudx(:,1:2,:)     = 1/dx * (-25/12*u(:,1:2,:) + 4*u(:,2:3,:) - 3*u(:,3:4,:) + 4/3*u(:,4:5,:) - 1/4*u(:,5:6,:));
        dudx(:,nx-1:nx,:) = -1/dx * (-25/12*u(:,nx-1:nx,:) + 4*u(:,nx-2:nx-1,:) - 3*u(:,nx-3:nx-2,:) + 4/3*u(:,nx-4:nx-3,:) - 1/4*u(:,nx-5:nx-4,:));
        
        % 4th order Central difference on 3rd cells
        
        dudx(:,3:nx-2,:) = 1/dx * ( rkw(2,3) * u(:,1:nx-4,:) + rkw(2,4) * u(:,2:nx-3,:) ...
            + rkw(2,5) * u(:,3:nx-2,:) + rkw(2,6) * u(:,4:nx-1,:) + rkw(2,7) * u(:,5:nx,:) );
        
    elseif order == 6
        
        % 6th order Upwind on boundaries
        
        dudx(:,1:3,:)     = +1/dx * (upw(6,1)*u(:,1:3,:) + upw(6,2)*u(:,2:4,:) + upw(6,3)*u(:,3:5,:) + upw(6,4)*u(:,4:6,:) + upw(6,5)*u(:,5:7,:) + upw(6,6)*u(:,6:8,:) + upw(6,7)*u(:,7:9,:));
        dudx(:,nx-2:nx,:) = -1/dx * (upw(6,1)*u(:,nx-2:nx,:) + upw(6,2)*u(:,nx-3:nx-1,:) + upw(6,3)*u(:,nx-4:nx-2,:) + upw(6,4)*u(:,nx-5:nx-3,:) + upw(6,4)*u(:,nx-5:nx-3,:) + upw(6,5)*u(:,nx-6:nx-4,:) + upw(6,6)*u(:,nx-7:nx-5,:) + upw(6,7)*u(:,nx-8:nx-6,:));
        
        % 2nd order Central difference on 2nd cells
        %     dudx(:,[2 nx-1],:) = 1/dx * ( rkw(1,4) * u(:,[1 nx-2],:) + rkw(1,6) * u(:,[3 nx],:) );
        
        % 4th order Central difference on 3rd cells
        
        dudx(:,3:nx-2,:) = 1/dx * ( rkw(2,3) * u(:,1:nx-4,:) + rkw(2,4) * u(:,2:nx-3,:) ...
            + rkw(2,5) * u(:,3:nx-2,:) + rkw(2,6) * u(:,4:nx-1,:) + rkw(2,7) * u(:,5:nx,:) );
        
    elseif order == 400
        
        % 2nd order Upwind on boundaries
        
        dudx(:,1,:)  = +1/dx * (upw(2,1)*u(:,1,:) + upw(2,2)*u(:,2,:) + upw(2,3)*u(:,3,:));
        dudx(:,nx,:) = -1/dx * (upw(2,1)*u(:,nx,:) + upw(2,2)*u(:,nx-1,:) + upw(2,3)*u(:,nx-2,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudx(:,[2 nx-1],:) = 1/dx * ( rkw(1,4) * u(:,[1 nx-2],:) + rkw(1,6) * u(:,[3 nx],:) );
        
        % 4th order Central difference on 3rd cells
        
        dudx(:,3:nx-2,:) = 1/dx * ( rkw(2,3) * u(:,1:nx-4,:) + rkw(2,4) * u(:,2:nx-3,:) ...
            + rkw(2,5) * u(:,3:nx-2,:) + rkw(2,6) * u(:,4:nx-1,:) + rkw(2,7) * u(:,5:nx,:) );
        
    elseif order == 4000
        
        % 1st order Upwind on boundaries
        
        dudx(:,1,:)  = 1/1/dx * (u(:,2,:) - u(:,1,:));
        dudx(:,nx,:) = 1/1/dx * (u(:,nx,:) - u(:,nx-1,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudx(:,[2 nx-1],:) = 1/dx * ( rkw(1,4) * u(:,[1 nx-2],:) + rkw(1,6) * u(:,[3 nx],:) );
        
        % 4th order Central difference on 3rd cells
        
        dudx(:,3:nx-2,:) = 1/dx * ( rkw(2,3) * u(:,1:nx-4,:) + rkw(2,4) * u(:,2:nx-3,:) ...
            + rkw(2,5) * u(:,3:nx-2,:) + rkw(2,6) * u(:,4:nx-1,:) + rkw(2,7) * u(:,5:nx,:) );
        
    elseif order == 888 % richardson
        
        % 1st order Upwind on boundaries
        
        dudx(:,1,:)  = 1/1/dx * (u(:,2,:) - u(:,1,:));
        dudx(:,nx,:) = 1/1/dx * (u(:,nx,:) - u(:,nx-1,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudx(:,2:nx-1,:) = 1/dx * ( rkw(1,4) * u(:,1:nx-2,:) + rkw(1,6) * u(:,3:nx,:) );
        
        % Least squares
        
        dudx(:,3:nx-2,:) = 1/12/dx * ( ...
            1 * u(:,1:nx-4,:) + ...
            -8 * u(:,2:nx-3,:) + ...
            0 * u(:,3:nx-2,:) + ...
            8 * u(:,4:nx-1,:) + ...
            -1 * u(:,5:nx,:) );
        
    elseif order == 999 % least squares
        
        % 1st order Upwind on boundaries
        
        dudx(:,1,:)  = 1/1/dx * (u(:,2,:) - u(:,1,:));
        dudx(:,nx,:) = 1/1/dx * (u(:,nx,:) - u(:,nx-1,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudx(:,2:nx-1,:) = 1/dx * ( rkw(1,4) * u(:,1:nx-2,:) + rkw(1,6) * u(:,3:nx,:) );
        
        % Least squares
        
        dudx(:,3:nx-2,:) = 1/10/dx * ( ...
            -2 * u(:,1:nx-4,:) + ...
            -1 * u(:,2:nx-3,:) + ...
            0 * u(:,3:nx-2,:) + ...
            1 * u(:,4:nx-1,:) + ...
            2 * u(:,5:nx,:) );
        
    elseif order == 60
        
        % Upwind on boundaries
        
        dudx(:,1,:)  = 1/1/dx * (u(:,2,:) - u(:,1,:));
        dudx(:,nx,:) = 1/1/dx * (u(:,nx,:) - u(:,nx-1,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudx(:,[2 nx-1],:) = 1/dx * ( rkw(1,4) * u(:,[1 nx-2],:) + rkw(1,6) * u(:,[3 nx],:) );
        
        % 4th order Central difference on 3rd cells
        
        dudx(:,[3 nx-2],:) = 1/dx * ( rkw(2,3) * u(:,[1 nx-4],:) + rkw(2,4) * u(:,[2 nx-3],:) ...
            + rkw(2,5) * u(:,[3 nx-2],:) + rkw(2,6) * u(:,[4 nx-1],:) + rkw(2,7) * u(:,[5 nx],:) );
        
        % 6th order Central difference on interior cells
        
        dudx(:,4:nx-3,:) = 1/dx * ( rkw(3,2) * u(:,1:nx-6,:) + rkw(3,3) * u(:,2:nx-5,:) + rkw(3,4) * u(:,3:nx-4,:) ...
            + rkw(3,5) * u(:,4:nx-3,:) + rkw(3,6) * u(:,5:nx-2,:) + rkw(3,7) * u(:,6:nx-1,:) + rkw(3,8) * u(:,7:nx,:) );
        
    elseif order == 80
        
        % Upwind on boundaries
        
        dudx(:,1,:)  = 1/1/dx * (u(:,2,:) - u(:,1,:));
        dudx(:,nx,:) = 1/1/dx * (u(:,nx,:) - u(:,nx-1,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudx(:,[2 nx-1],:) = 1/dx * ( rkw(1,4) * u(:,[1 nx-2],:) + rkw(1,6) * u(:,[3 nx],:) );
        
        % 4th order Central difference on 3rd cells
        
        dudx(:,[3 nx-2],:) = 1/dx * ( rkw(2,3) * u(:,[1 nx-4],:) + rkw(2,4) * u(:,[2 nx-3],:) ...
            + rkw(2,5) * u(:,[3 nx-2],:) + rkw(2,6) * u(:,[4 nx-1],:) + rkw(2,7) * u(:,[5 nx],:) );
        
        % 6th order Central difference on interior cells
        
        dudx(:,[4 nx-3],:) = 1/dx * ( rkw(3,2) * u(:,[1 nx-6],:) + rkw(3,3) * u(:,[2 nx-5],:) + rkw(3,4) * u(:,[3 nx-4],:) ...
            + rkw(3,5) * u(:,[4 nx-3],:) + rkw(3,6) * u(:,[5 nx-2],:) + rkw(3,7) * u(:,[6 nx-1],:) + rkw(3,8) * u(:,[7 nx],:) );
        
        % 8th order Central difference on interior cells
        
        dudx(:,5:nx-4,:) = 1/dx * ( ...
            rkw(4,1) * u(:,1:nx-8,:) + ...
            rkw(4,2) * u(:,2:nx-7,:) + ...
            rkw(4,3) * u(:,3:nx-6,:) + ...
            rkw(4,4) * u(:,4:nx-5,:) + ...
            rkw(4,5) * u(:,5:nx-4,:) + ...
            rkw(4,6) * u(:,6:nx-3,:) + ...
            rkw(4,7) * u(:,7:nx-2,:) + ...
            rkw(4,8) * u(:,8:nx-1,:) + ...
            rkw(4,9) * u(:,9:nx,:) );
        
    elseif order == 15
        
        % 1st order Upwind on boundaries
        
        dudx(:,1,:)  = 1/1/dx * (u(:,2,:) - u(:,1,:));
        dudx(:,nx,:) = 1/1/dx * (u(:,nx,:) - u(:,nx-1,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudx(:,[2 nx-1],:) = 1/dx * ( rkw(1,4) * u(:,[1 nx-2],:) + rkw(1,6) * u(:,[3 nx],:) );
        
        % 4th order Central difference on 3rd cells
        
        dudx(:,3:nx-2,:) = 1/dx/10 * ( ...
            -2 * u(:,1:nx-4,:) + ...
            -1 * u(:,2:nx-3,:) + ...
            +1 * u(:,4:nx-1,:) + ...
            +2 * u(:,5:nx,:) );
        
    elseif order == 123
        
        % 1st order Upwind on boundaries
        
        dudx(:,1,:)  = 1/1/dx * (u(:,2,:) - u(:,1,:));
        dudx(:,nx,:) = 1/1/dx * (u(:,nx,:) - u(:,nx-1,:));
        
        % 2nd order Central difference on 2nd cells
        
        dudx(:,2:nx-1,:) = 1/dx * ( rkw(1,4) * u(:,1:nx-2,:) + rkw(1,6) * u(:,3:nx,:) );
        
        % CR4s on interior cells
        
        A  = 1239;
        A1 = 272;
        A2 = 1036;
        A4 = 0;
        A8 = -69;
        
        dudx(:,9:nx-8,:) = 1/A * ( ...
            A1/2/1/dx * (u(:,9+1:nx-8+1,:) - u(:,9-1:nx-8-1,:)) + ...
            A2/2/2/dx * (u(:,9+2:nx-8+2,:) - u(:,9-2:nx-8-2,:)) + ...
            A4/2/4/dx * (u(:,9+4:nx-8+4,:) - u(:,9-4:nx-8-4,:)) + ...
            A8/2/8/dx * (u(:,9+8:nx-8+8,:) - u(:,9-8:nx-8-8,:)) ...
            );
        
    else
        
        error('select derivative scheme')
        
    end
    
    
end