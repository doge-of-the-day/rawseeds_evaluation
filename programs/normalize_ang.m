function angles = normalize_ang(angles)
%--------------------------------------------------------------------------
%                  RAWSEEDS METRICS COMPUTATION TOOLKIT
%                       http://www.rawseeds.org/
%--------------------------------------------------------------------------
%   function   angles = normalize_ang(angles
%
%   Normalizes angles to interval (-pi,pi]
%
%   Authors: J. Neira, J. Tardos
%--------------------------------------------------------------------------
%   Version: 1.0    jul-2002
%--------------------------------------------------------------------------
%   History:
%--------------------------------------------------------------------------

angles = angles + (2*pi)*floor((pi-angles)/(2*pi));
