function tac=trel(t1,t2)
%--------------------------------------------------------------------------
%                  RAWSEEDS METRICS COMPUTATION TOOLKIT
%--------------------------------------------------------------------------
%   function   tac=trel(t1,t2)
%
%   Computes the relative transformations: inv(t1)*t2
%
%   Authors: J. Neira, J. Tardos
%--------------------------------------------------------------------------
%   Version: 1.1    oct-2009
%--------------------------------------------------------------------------
%   History: 
%            1.0    jul-2002 Original
%            1.1    oct-2009 vectorized to gain efficiency  J.D.Tardos
%--------------------------------------------------------------------------

if size(t1,1) ~= 3,
   error('TREL: t1 is not a transformation!!!');
end;

if size(t2,1) ~= 3,
   error('TREL: t2 is not a transformation!!!');
end;

s = sin(t1(3,:));
c = cos(t1(3,:));
dx = t2(1,:)-t1(1,:);
dy = t2(2,:)-t1(2,:);
tac =  [ c .* dx + s .* dy
        -s .* dx + c .* dy
        normalize_ang(t2(3,:)-t1(3,:)) ];
         
