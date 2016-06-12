function [t_rpe,r_rpe]=get_RPE(X,GT_rel,plots)
%--------------------------------------------------------------------------
%                  RAWSEEDS METRICS COMPUTATION TOOLKIT
%                       http://www.rawseeds.org/
%--------------------------------------------------------------------------
%   function [T-RPE,R-RPE] = get_RPE(X,GT_rel)
% 
%   Computes the Relative Pose Error:
%     x_ij = inv(x_i)    * x_j
%     e_ij = inv(xGT_ij) * x_ij
%     t_rpe = sqrt(sum(||trans(e_ij)||^2)/n);
%     r_rpe = sqrt(sum(rot(e_ij)^2)/n)*180/pi;
%     
%   Authors: C. Cadena & J.D.Tardos, University of Zaragoza, Spain
%--------------------------------------------------------------------------
%   Version: 1.0    oct-2009
%--------------------------------------------------------------------------
%   History:
%-------------------------------------------------------------------------- 
% 
% inputs:  
%   X = [t,x,y,theta], matrix n x 4, n steps, from the SLAM solution,
%          t is the vector of timestamps in each step.
%
%   GT_rel = [t_i, t_j, dx, dy, dtheta], matrix m x 5, with the relations
%          between the timestamp t_i and the timestamp t_j.
%   
% outputs: 
%   T-RPE = translational relative pose error in meters
%   R-RPE = rotational relative pose error in degrees
%
% [...] = get_RPE(X,,GT_rel,plots), , plots = 1 to plot results
%--------------------------------------------------------------------------

%% global parameters
global PARAMETERS;

%% checking inputs
if nargin < 2
    error('Not enough input arguments.'); 
elseif size(X,2)<4
    error('X must be of size n x 4.');
end
if nargin <3
    plots=0;
end

%% Remove GT values before and after the SLAM timestamps
tini = X(1,1);
tend = X(end,1);

t_i=GT_rel(:,1);
t_j=GT_rel(:,2);
xGT_ij=GT_rel(:,3:5);

good= (t_i>=tini) & (t_i<=tend) & (t_j>=tini) & (t_j<=tend);
t_i = t_i(good);
t_j = t_j(good);
xGT_ij = xGT_ij(good,:);

n=size(xGT_ij,1);

%Interpolate SLAM trajectory at the GT timestamps
X(:,4)=unwrap(X(:,4));
x_i = interp1(X(:,1),X(:,2:4),t_i);
x_j = interp1(X(:,1),X(:,2:4),t_j);

% Compute Relative Pose Error
x_ij = trel(x_i', x_j')';
e_ij = trel(xGT_ij', x_ij')';

tvector = sum(e_ij(:,1:2).^2,2);
rvector = e_ij(:,3).^2;
t_rpe = sqrt(mean(tvector));
r_rpe = sqrt(mean(rvector))*180/pi;

%% plots
if plots
    DATASET=PARAMETERS.DatasetName;    
    SOLname=PARAMETERS.SolutionName(1:end-4);
        
    figure('name','RAWSEEDS METRICS COMPUTATION TOOLKIT')
    subplot(2,1,1), plot(sqrt(tvector)); ylabel('T-RPE [m]'); 
    xlabel('relation number')    
    title([DATASET,': RPE with ',SOLname],'Interpreter','none')
    subplot(2,1,2), plot(sqrt(rvector)*180/pi); ylabel('R-RPE [deg]'); 
    xlabel('relation number')    
    
    [st_i, i_i]=sort(t_i);
    figure('name','RAWSEEDS METRICS COMPUTATION TOOLKIT')
    subplot(2,1,1), plot(st_i-tini,sqrt(tvector(i_i))); ylabel('T-RPE [m]'); 
    xlabel('t_i')
    title([DATASET,': RPE with ',SOLname,' sorted by t_i'],'Interpreter','none')
    subplot(2,1,2), plot(st_i-tini,sqrt(rvector(i_i))*180/pi); ylabel('R-RPE [deg]'); 
    xlabel('t_i')
    
    [st_j,i_j]=sort(t_j);
    figure('name','RAWSEEDS METRICS COMPUTATION TOOLKIT')    
    subplot(2,1,1), plot(st_j-tini,sqrt(tvector(i_j))); ylabel('T-RPE [m]'); 
    xlabel('t_j')        
    title([DATASET,': RPE with ',SOLname,' sorted by t_j'],'Interpreter','none')
    subplot(2,1,2), plot(st_j-tini(1),sqrt(rvector(i_j))*180/pi); ylabel('R-RPE [deg]'); 
    xlabel('t_j')
    
end
