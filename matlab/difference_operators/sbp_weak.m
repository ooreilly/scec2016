function [xp,xm,Pp,Pm,Qp,Qm] = sbp_weak(order,n,h,test)
% [xp,xm,Pp,Pm,Qp,Qm] = sbp_weak(order,n,h)
%
% Construct SBP staggered grid operators
% that satisfy the property (Q+)_00  = (Q_-)_00,  (H+)_00 = (H-)_00
%
% Input:
% order : Order of accuracy
% n     : Controls number of grid points
% h     : Grid spacing
%
% Output
% xp,xm         : Grid vectors xp (n+1 grid points) xm (n+2 grid points)
%                 xp = [0 h 2h + ... ], xm = [0 h/2 3h/2 +  ... ]
% Pp,Pm,Qp,Qm,  : Staggered grid operators

if nargin < 4
  test = false;
end

    switch order
      case 4
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_weak_4th(n,h,test);
      case 6
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_weak_6th(n,h,test);
      otherwise
       error('SBP staggered grid operator not implemented');   
      end

xp = xp';
xm = xm';


end
