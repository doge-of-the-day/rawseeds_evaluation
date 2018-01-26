%--------------------------------------------------------------------------
%                 RAWSEEDS METRICS COMPUTATION TOOLKIT
%                       http://www.rawseeds.org/
%--------------------------------------------------------------------------
%   benchmark.m
%
%   This script loads a SLAM benchmark solution and the ground-truth files
%   from the selected RAWSEEDS dataset and computes two metrics:
%   - Absolute Trajectory Error (ATE)
%   - Relative Pose Error (RPE) 
%   For more information, refer to the project site: www.rawseeds.org
%
%   Authors: C. Cadena & J.D.Tardos, University of Zaragoza, Spain
%--------------------------------------------------------------------------
%   Version:    1.0    31-oct-2009
%               1.1    17-dec-2009
%--------------------------------------------------------------------------
%   History:    1.1    added choice of interval used for GT alignment in ATE
%--------------------------------------------------------------------------

close all
clear all

global PARAMETERS

PARAMETERS.DataPath      = '../datasets';

PARAMETERS.ComputeATE    = 1;
PARAMETERS.ComputeRPE    = 1;
PARAMETERS.show_plots    = 1;
PARAMETERS.scale_unknown = 0;
PARAMETERS.alignment_path= NaN;  % interval used for GT alignment in ATE
                                 %  - [1:100] to use first 100 steps
                                 %  - NaN to use full trajectory

% selection of the dataset to use
% PARAMETERS.DatasetName   = 'Bicocca_2009-02-25a';  % Indoor
% PARAMETERS.DatasetName   = 'Bicocca_2009-02-25b';  % Indoor
 PARAMETERS.DatasetName   = 'Bicocca_2009-02-26a';  % Indoor
% PARAMETERS.DatasetName   = 'Bicocca_2009-02-26b';  % Indoor
% PARAMETERS.DatasetName   = 'Bicocca_2009-02-27a';  % Indoor
% PARAMETERS.DatasetName   = 'Bovisa_2008-09-01';    % Mixed
% PARAMETERS.DatasetName   = 'Bovisa_2008-10-06';    % Mixed
% PARAMETERS.DatasetName   = 'Bovisa_2008-10-11a';   % Mixed
% PARAMETERS.DatasetName   = 'Bovisa_2008-10-04';    % Outdoor
% PARAMETERS.DatasetName   = 'Bovisa_2008-10-07';    % Outdoor
% PARAMETERS.DatasetName   = 'Bovisa_2008-10-11b';   % Outdoor

% selection of the algorithm (Benchmark Solution) to be assessed
% PARAMETERS.SolutionName  = 'CI-GRAPH Stereo.csv';   % Bicocca_2009-02-25b
% PARAMETERS.SolutionName  = 'GraphSLAM Laser.csv';   % Bicocca_2009-02-25b
% PARAMETERS.SolutionName  = 'H-SLAM Trinocular.csv'; % Bicocca_2009-02-25b
% PARAMETERS.SolutionName  = 'EKF Monocular.csv';     % Bovisa_2008-10-04
%PARAMETERS.SolutionName  = 'Bicocca_2009-02-25b-GROUNDTRUTH_CARTOGRAPHER.csv';
PARAMETERS.SolutionName  = 'muse_beam_model_default.csv';


% selection of the ground truth to be used as a reference
% PARAMETERS.GT_name       = '-GROUNDTRUTH.csv';  % only for Bicocca/indoor
%  PARAMETERS.GT_name       = '-GT-extended.csv';  % only for Bicocca/indoor
  PARAMETERS.GT_name       = '-GROUNDTRUTH-CARTOGRAPHER.csv';  % only for Bicocca/indoor
% PARAMETERS.GT_name       = '-GPS.csv';          % only for Bovisa/outdoor
  PARAMETERS.GT_rel_name   = '-GROUNDTRUTH-CARTOGRAPHER.relations';        % for all datasets

example = 0;  % Change this number to try the different examples provided, 
              % instead of using the previous parameters 
switch example
  case 1      
    PARAMETERS.DatasetName   = 'Bicocca_2009-02-25b';  % Indoor
    PARAMETERS.SolutionName  = 'GraphSLAM Laser.csv';
    PARAMETERS.GT_name       = '-GT-extended.csv';
  case 2      
    PARAMETERS.DatasetName   = 'Bicocca_2009-02-25b';  % Indoor
    PARAMETERS.SolutionName  = 'CI-GRAPH Stereo.csv';
    PARAMETERS.GT_name       = '-GT-extended.csv';
  case 3      
    PARAMETERS.DatasetName   = 'Bicocca_2009-02-25b';  % Indoor
    PARAMETERS.SolutionName  = 'H-SLAM Trinocular.csv';
    PARAMETERS.GT_name       = '-GT-extended.csv';
  case 4
    PARAMETERS.DatasetName   = 'Bovisa_2008-10-04';    % Outdoor
    PARAMETERS.SolutionName  = 'EKF Monocular.csv';
    PARAMETERS.GT_name       = '-GPS.csv';
end

%% import file of solution
folder=[PARAMETERS.DataPath,'/',PARAMETERS.DatasetName,'/'];
dataset=PARAMETERS.DatasetName;

disp('--------------------------------------------------------------------');
disp('               RAWSEEDS METRICS COMPUTATION TOOLKIT                 ');
disp('--------------------------------------------------------------------');
fprintf('Loading solution: %s \n\n',[folder,PARAMETERS.SolutionName]);
X=importdata([folder,PARAMETERS.SolutionName]);

if PARAMETERS.ComputeATE
    disp('--------------------------------------------------------------------');
    disp('                    Absolute Trajectory Error                       ');
    disp('--------------------------------------------------------------------');
    % import file of GT
    if strcmpi(PARAMETERS.GT_name(2:4),'GPS')
%        file_GT=[dataset(1:find(dataset=='_',1,'last')-1),PARAMETERS.GT_name];
        file_GT=[dataset,PARAMETERS.GT_name];
        fprintf('Loading GT file:  %s \n',[folder,file_GT]);
        GT=read_gps([folder,file_GT],[X(1,1),X(end,1)]);        
    else
        file_GT=[dataset,PARAMETERS.GT_name];
        fprintf('Loading GT file:  %s \n',[folder,file_GT]);
        GT=importdata([folder,file_GT]);
        GT=GT(:,1:4);
    end
    
    % computing ATE
    ATE=get_ATE(X,GT,PARAMETERS.show_plots,PARAMETERS.scale_unknown);    
    fprintf('Dataset:      %s \n',PARAMETERS.DatasetName);
    fprintf('Solution:     %s \n',PARAMETERS.SolutionName(1:end-4));
    fprintf('Ground Truth: %s \n',PARAMETERS.GT_name(2:end-4));
    fprintf('ATE (mean)      : %8.4f m \n',ATE(1));
    fprintf('ATE (std)       : %8.4f m \n',ATE(2));
    fprintf('ATE (mean-3*std): %8.4f m \n',ATE(3));
    fprintf('ATE (mean+3*std): %8.4f m \n\n',ATE(4));
end

if PARAMETERS.ComputeRPE
    disp('--------------------------------------------------------------------');
    disp('                      Relative Pose Error                           ');
    disp('--------------------------------------------------------------------');
    % import file of relations
    file_rel=[dataset,PARAMETERS.GT_rel_name];
    fprintf('Loading GT file:  %s \n',[folder,file_rel]);
    GT_rel=importdata([folder,file_rel]);
    GT_rel=GT_rel(:,[1:4,8]);
    
    % computing RPE
    [t_rpe,r_rpe]=get_RPE(X,GT_rel,PARAMETERS.show_plots);
    fprintf('Dataset:  %s \n',PARAMETERS.DatasetName);
    fprintf('Solution: %s \n',PARAMETERS.SolutionName(1:end-4));
    fprintf('T-RPE: %8.4f m \nR-RPE: %8.4f deg \n\n',t_rpe,r_rpe);
end
