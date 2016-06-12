function tac=tcomp(tab,tbc)
%--------------------------------------------------------------------------
%                  RAWSEEDS METRICS COMPUTATION TOOLKIT
%--------------------------------------------------------------------------
%   function   tac=tcomp(tab,tbc)
%
%   composes two transformations in 2D
%
%   Authors: J. Neira, J. Tardos
%--------------------------------------------------------------------------
%   Version: 1.2    oct-2009
%--------------------------------------------------------------------------
%   History: 
%            1.0    jul-2002 Original
%            1.1    jan-2008 tbc can be 3xn now, n > 1, C. Cadena
%            1.2    oct-2009 vectorized to gain efficiency  J.D.Tardos
%--------------------------------------------------------------------------

if size(tab,1) ~= 3,
   error('TCOMP: tab is not a transformation!!!');
end;

if size(tbc,1) ~= 3,
   error('TCOMP: tbc is not a transformation!!!');
end;

s = sin(tab(3));
c = cos(tab(3));

tac =  [tab(1) + [c -s]*tbc(1:2,:)
        tab(2) + [s  c]*tbc(1:2,:)
        normalize_ang(tab(3)+tbc(3,:)) ];
         
