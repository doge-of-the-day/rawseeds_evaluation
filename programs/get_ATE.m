function [ATE,errors,X_e,GT_e,RT]=get_ATE(X,GT,plots,unknown_scale)
%--------------------------------------------------------------------------
%                  RAWSEEDS METRICS COMPUTATION TOOLKIT
%                       http://www.rawseeds.org/
%--------------------------------------------------------------------------
%   function [ATE,e,RT] = get_ATE(X,GT)
%
%   Aligns the trajectory X with the ground-truth trajectory GT and
%   computes the Absolute Trajectory Error:
%     e_j = inv(xGT_j) * x_j
%     d_j = ||trans(e_j)||
%     d = mean(d_j)
%     s = std_dev(d_j)
%     ATE = [d, s, d-3*s, d+3*s]
%
%   Authors: C. Cadena & J.D.Tardos, University of Zaragoza, Spain
%
%--------------------------------------------------------------------------
%   Version: 1.0    oct-2009
%            1.1    nov-2009
%--------------------------------------------------------------------------
%   History: 1.1 added alignment_path parameters
%--------------------------------------------------------------------------
%
% inputs:
%   X = [t,x,y,theta], matrix n x 4, n steps, from the SLAM solution,
%          t is the vector of timestamps in each step.
%   GT= [tGT,xGT,yGT,thetaGT], matrix m x 4, m steps, from the Ground Truth,
%          tGT is the vector of timestamps in each step.
%   plots = 1 if you want to plot errors
%   unknows_scale = 1 if your SLAm algorithm does not recover the scale
%        (for example, it uses pure monocular SLAM, without odometry).
%        In this case, this function computes the transformation AND the
%        scale that better aligns your solution with the ground truth.
%
% outputs:
%   ATE =[d,s,a,b],
%       d = mean of translation error
%       s = standard deciation of the translation errors
%       a,b = extremes of the 3s confidence interval of the translation
%       errors
%   e = matrix n x 3 with the errors in x,y and theta for indoors datasets
%     = matrix n x 2 with the errors in x and y for mixed and outdoors datasets
%   RT= Vector with the transformation from your solution to the GT
%
%--------------------------------------------------------------------------

global PARAMETERS;

%% Remove GT values before and after the SLAM timestamps
tini=X(1,1);
tend=X(end,1);
good= (GT(:,1)>=tini) & (GT(:,1)<=tend);
GT = GT(good,:);
t_GT = GT(:,1)-tini;

%Interpolate SLAM trajectory at the GT timestamps
X(:,4)=unwrap(X(:,4));
X_e = interp1(X(:,1),X(:,2:4),GT(:,1));
GT_e= GT(:,2:end);

%% optimize alignment
disp('Aligning trajectory with GT...');
optioptions = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
                  'Display','off');
if unknown_scale
    loc_ini=[0 0 0 1];
    if(isnan(PARAMETERS.alignment_path))
        loc_f=lsqnonlin('scale_alignment_error',loc_ini,[-inf -inf -pi 0],[Inf Inf pi Inf],...
            optioptions,...
            GT_e(:,1:2)',X_e(:,1:2)');
    else
        loc_f=lsqnonlin('scale_alignment_error',loc_ini,[-inf -inf -pi 0],[Inf Inf pi Inf],...
            optioptions,...
            GT_e(PARAMETERS.alignment_path,1:2)',X_e(PARAMETERS.alignment_path,1:2)');
    end        
    s=loc_f(4);
    X(:,2:4)=tcomp(loc_f(1:3)',[s*X(:,2:3)   X(:,4)]')';
    X_e     =tcomp(loc_f(1:3)',[s*X_e(:,1:2) X_e(:,3)]')';
else
    loc_ini=[0 0 0];
    if(isnan(PARAMETERS.alignment_path))
        loc_f=lsqnonlin('alignment_error',loc_ini,[-inf -inf -pi],[Inf Inf pi],...
            optioptions,...
            GT_e(:,1:2)',X_e(:,1:2)');    
    else
        loc_f=lsqnonlin('alignment_error',loc_ini,[-inf -inf -pi],[Inf Inf pi],...
            optioptions,...
            GT_e(PARAMETERS.alignment_path,1:2)',X_e(PARAMETERS.alignment_path,1:2)');    
    end
    
    X(:,2:4)=tcomp(loc_f',X(:,2:4)')';
    X_e=tcomp(loc_f',X_e')';
end
RT=loc_f(1:3)';
  fprintf('Aligning transformation: [%3.4f %3.4f %3.4f] \n',RT);
if unknown_scale
  fprintf('Aligning scale factor:   %3.4f \n',s);
end
fprintf('\n');

%% Compute error
n = size(GT_e,2);
errors=X_e(:,1:n)-GT_e(:,1:n);
if n>2, errors(:,3)=normalize_ang(errors(:,3)); end

error_modulus=sqrt(sum(errors(:,1:2)'.^2))';
if n>2, error_modulus(:,2)=abs(errors(:,3)')'; end

e_mean=mean(error_modulus);
e_std=std(error_modulus);
e_a=e_mean-3*e_std;
e_b=e_mean+3*e_std;

% e_max=max(error_modulus);
% fprintf('\nMean of Absolute Error:\n%2.5f m \t',e_mean(1));
% if n>2, fprintf('%2.5f rad',e_mean(2)); end
% fprintf('\nMaximum Absolute Error:\n%2.5f m \t',e_max(1));
% if n>2, fprintf('%2.5f rad \n',e_max(2)); end
% fprintf(' \n');


%% Form ATE vector
md=mean(error_modulus(:,1));
sd=std(error_modulus(:,1));
ATE =[md, sd, md-3*sd, md+3*sd];

%% plots
if plots
    DATASET=PARAMETERS.DatasetName;
    GTname=PARAMETERS.GT_name(2:end-4);
    SOLname=PARAMETERS.SolutionName(1:end-4);
    figure('name','RAWSEEDS METRICS COMPUTATION TOOLKIT')
    subplot(n,1,1)
    plot(X(:,1)-tini,X(:,2),'b.',t_GT,GT(:,2),'r.')
    title(DATASET,'Interpreter','none');
    legend(SOLname,GTname)
    ylabel('x [m]')
    subplot(n,1,2)
    plot(X(:,1)-tini,X(:,3),'b.',t_GT,GT(:,3),'r.')
    ylabel('y [m]')
    if n>2
        subplot(n,1,n)
        plot(X(:,1)-tini,X(:,4),'b.',t_GT,GT(:,n+1),'r.')
        ylabel('\theta [rad]')
    end
    xlabel('time [s]')

    figure('name','RAWSEEDS METRICS COMPUTATION Errors')
    subplot(n,1,1)
    plot(t_GT, errors(:,1),'b.')
    ylabel('error in x [m]')
    title([DATASET,': Errors ',SOLname,' vs ',GTname],'Interpreter','none');
    subplot(n,1,2)
    plot(t_GT, errors(:,2),'b.')
    ylabel('error in y [m]')
    if n>2
        subplot(n,1,n)
        plot(t_GT, errors(:,n),'b.')
        ylabel('error in \theta [rad]')
    end

    xlabel('time [s]')

    X_e(1,:)
    GT_e(1,:)
    figure('name','RAWSEEDS METRICS COMPUTATION TOOLKIT')
plot(X(:,2),X(:,3),'c',X_e(:,1),X_e(:,2),'bx',GT_e(:,1),GT_e(:,2),'-r.');
    hold on;
    xlabel('x [m]'), ylabel('y [m]'), axis equal
    title(DATASET,'Interpreter','none')
    legend(SOLname,[SOLname,' in the GT zone'],GTname)

    figure('name','Our histogramm')
    hist3(errors(:,1:2), [100,100]);
    pointsize = 10;
    scatter3(errors(:,1), errors(:,2),errors(:,3), pointsize)

    figure('name','error distribution')
    hist3(errors(:,1:2), [100,100]);


    figure('name','RAWSEEDS METRICS COMPUTATION TOOLKIT')
    [yh,xh]=hist(error_modulus(:,1),100);
    bar(xh,yh,'c','BarWidth',0.8,'EdgeColor','c');
    hold on;
    plot([e_mean(1), e_mean(1)],[0,max(yh)],'b',[e_a(1), e_a(1)],[0,max(yh)],'r',[e_b(1), e_b(1)],[0,max(yh)],'r');
    xlabel('Position Error [m]'), ylabel('Frequency')
    title([DATASET,': ATE, ',SOLname,' vs ',GTname],'Interpreter','none')
    text(e_mean(1),max(yh),' Mean','Color','b'), text(e_a(1),max(yh)/2,' -3\sigma','Color','r'),text(e_b(1),max(yh)/2,' +3\sigma','Color','r');
end
