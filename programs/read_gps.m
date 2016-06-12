function GT = read_gps(file,t)
%--------------------------------------------------------------------------
%                  RAWSEEDS METRICS COMPUTATION TOOLKIT
%                       http://www.rawseeds.org/
%--------------------------------------------------------------------------
%   function   GT = read_gps(file,t)
%
%   Imports a GPS file and transforms the geographical coordinates to
%   meters in 2D
%
%   Authors: J. Civera, C. Cadena, University of Zaragoza, Spain
%--------------------------------------------------------------------------
%   Version: 1.0    oct-2009
%--------------------------------------------------------------------------
%   History:
%--------------------------------------------------------------------------

gps_data = importdata(file);

useful_gps_data = [];
for i=1:size(gps_data,1)
    timestamp = str2double(gps_data{i}(1:17));
    c=find(gps_data{i}==',');
    if gps_data{i}(25)=='A'
        lat = str2double(gps_data{i}(c(3)+1:c(4)-1));
        lon = str2double(gps_data{i}(c(5)+1:c(6)-1));
        height = str2double(gps_data{i}(c(10)+1:c(11)-1));
    elseif gps_data{i}(25)=='T'
        latitud_error = str2double(gps_data{i}(c(7)+1:c(8)-1));
        longitude_error = str2double(gps_data{i}(c(8)+1:c(9)-1));
        if timestamp>t(1)&& timestamp<t(2)
            if (latitud_error<0.2)&&(longitude_error<0.2)
                useful_gps_data = [useful_gps_data; [timestamp, ...
                     floor(lat/100)+mod(lat/100,floor(lat/100))*100/60, ...
                     floor(lon/100)+mod(lon/100,floor(lon/100))*100/60],...
                     height];        
            end
        end
    end    
end
ellipsoid = w_ellipsoid(6378137, 6356752.3142);
utm = geo2UTM(useful_gps_data(:,3), useful_gps_data(:,2), ellipsoid);

GT = [useful_gps_data(:,1), utm(:,1)-utm(1,1), utm(:,2)-utm(1,2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  ellipsoid = w_ellipsoid(s_major, s_minor)
%   function   ellipsoid = w_ellipsoid(s_major, s_minor)
%
%   Generates an ellipsoid since semi-principal axises
%
%   Author:  O. Grasa

% WGS84: semi-major axis = 6378137 
%        semi-minor axis = 6356752.3142 
a = s_major;
b = s_minor;
ep = sqrt( (a * a) - (b * b)) / b;
ep2 = ep * ep;
c = (a * a) / b;
ellipsoid = [a b ep ep2 c];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function utm = geo2UTM(lon, lat, ellipsoid)
%   function   [gt,ws,plots]=ATE_options(varargin)
%
%   Tranforms geographical coordinates to UTM
%
%   Author:  O. Grasa
latRad = (lat * pi) / 180.0;
lonRad = (lon * pi) / 180.0;
tzone = floor(lon / 6 + 31);
merCentral = ((tzone * 6 - 183) * pi) / 180;
ang_dist = lonRad - merCentral;
[x,y,h] = cs(latRad, ang_dist, ellipsoid);
utm = [x y tzone merCentral];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,h] = cs(latRad, ang_dist, ellipsoid)
%   function   [x y h] = CS(latRad, ang_dist, ellipsoid)
%
%   Implementation of Coticchia-Surace equations for change of geographic
%   coordinates to UTM
%
%   Author:  O. Grasa
c = ellipsoid(5);
ep2 = ellipsoid(4);
A = cos(latRad) .* sin(ang_dist);
dseda = 0.5 .* log((1.0 + A) ./ (1.0 - A));
eta = atan( tan(latRad) ./ cos(ang_dist) ) - latRad;
aux = cos(latRad);
aux = aux .* aux;
ipsilon = c * 0.9996 ./ sqrt(1.0 + ep2 * aux);
sigma = 0.5 * ep2 .* dseda .* dseda .* aux;
A1 = sin(2.0 * latRad);
A2 = A1 .* aux;
J2 = latRad + A1 * 0.5;
J4 = 0.75 * J2 + 0.25 * A2;
J6 = (5.0 * J4 + A2 .* aux)/3;
alpha = 0.75 * ep2;
aux = alpha .* alpha;
beta = (5.0 * aux) / 3.0;
gamma = (35.0 * aux .* alpha) / 27.0;
B = 0.9996 * c *(latRad - alpha .* J2 + beta * J4 - gamma * J6);
x = dseda .* ipsilon .* ( 1.0 + sigma / 3.0) + 500000.0;
y = eta .* ipsilon .* ( 1.0 + sigma) + B;
h = 'N';
if(y < 0)
    y = y + 10000000.0; 
    h = 'S';
end
