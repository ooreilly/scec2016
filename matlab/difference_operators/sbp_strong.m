function [xp,xm,Pp,Pm,Qp,Qm] = sbp_strong(order,n,h)
% [xp,xm,Pp,Pm,Qp,Qm] = sbp_weak(order,n,h)
%
% Input:
% order : Order of accuracy
% n     : Controls number of grid points
% h     : Grid spacing
%
% Output
% xp,xm         : Grid vectors xp (n+1 grid points) xm (n grid points)
%                 xp = [0 h 2h + ... ], xm = [h/2 3h/2 5h/2 +  ... ]
% Pp,Pm,Qp,Qm,  : Staggered grid operators

    switch order
      case 4
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_strong_4th(n,h);
      case 6
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_strong_6th(n,h);
      otherwise
       error('SBP staggered grid operator not implemented');   
    end

    Pm(:,1) = [];
    Pm(:,end) = [];
    Pm(1,:) = [];
    Pm(end,:) = [];
    Qp(:,1) = [];
    Qp(:,end) = [];
    Qm(1,:) = [];
    Qm(end,:) = [];   
    xm = xm(2:end-1);

    xp = xp';
    xm = xm';

end
