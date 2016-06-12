function d=alignment_error(vloc,p1,p2)
%--------------------------------------------------------------------------
%                  RAWSEEDS METRICS COMPUTATION TOOLKIT
%--------------------------------------------------------------------------
%   function   d = alignment_error(vloc,p1,p2)
%
%   Computes the euclidean distances between points p1 and vloc*p2, in 2D
%   vloc = (x, y, theta)' is a location vector
%
%   Authors: C. Cadena & J.D.Tardos, University of Zaragoza, Spain
%--------------------------------------------------------------------------
%   Version: 1.0    oct-2009
%--------------------------------------------------------------------------
%   History:
%--------------------------------------------------------------------------

s = sin(vloc(3));
c = cos(vloc(3));
p2aligned = [vloc(1) + [c -s]*p2
             vloc(2) + [s  c]*p2];
d=p1(:)-p2aligned(:);

