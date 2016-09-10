function [Dp Dm] = fs1(n,h)
% [Dp Dm] = fs1(n,h)
% Difference operators implementing the FS1 free-surface boundary condition 
% at each end-point.
% Dm is applied to the velocity and Dp is applied to the stress.
%
% Input:
%          n: number of grid points
%
% Output:
%    [Dp Dm]: Difference operator Dp of size (n + 1) x n and Dm of size n x (n + 1) 

% Construct Dp
n             = n + 1;
w             = 4; 
s             = rot90(vander(1:w))\((0:(w-1)).*(w/2-1/2+1).^([0 0:w-2]))';  
Dp            = spdiags(repmat(s(:)',[n+2 1]), -(w/2-1)-1:w/2-1, n+2, n+2);    
Dp            = Dp(2:end-1,2:end-2);
Dp(1,1)       = 2*s(3);
Dp(1,2)       = 2*s(4);
Dp(2,1)       = s(4) + s(2);
Dp(end,end)   = -Dp(1,1);
Dp(end,end-1) = -Dp(1,2);
Dp(end-1,end) = -Dp(2,1);
Dp            = Dp/h;

% Construct Dm
n             = n + 1;
w             = 4; 
s             = rot90(vander(1:w))\((0:(w-1)).*(w/2-1/2+1).^([0 0:w-2]))';  
Dm            = spdiags(repmat(s(:)',[n+2 1]), -(w/2-1)-1:w/2-1, n+2, n+2);    
Dm            = Dm(3:end-2,2:end-2);
Dm(1,2)       = s(3) + s(1);
Dm(end,end-1) = -Dm(1,2);
Dm            = Dm/h;
