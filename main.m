%% TODO:: Plotting / Subplots / Convergence
%% TODO:: DetermineLambdas / Datasets / HyperParameters
%%
clc;
clear;
close all;
format long;

%%

root_path = './';           % Root path
image_folder = 'plots/';    % Image folder
data_folder = 'datasets/';  % Dataset folder
backup_folder = 'backup/';  % Backup folder for auxiliary matrices
showfig = true;             % Show figure
% print_output = false;       % Print output
print_output = true;

%% 

d = 2;                      % Dimension of data points
p = [2,3,4,10];             % Num of clusters
N = [200,600,400,3000];     % Num of observations

handle_outlier = false;     % Outlier handling

%% Parameters

i = 2;

lambda = 0.05;
delta = 10^(-4);
L = 1 + N(i)*lambda/delta;
x0 = sparse(N(i),d);
% tol = 10^(-2);
tol = 10^(-1);
mu = 1;
epsilon = 10^(-2);

load([root_path,data_folder,'data_',num2str(i),'_size',num2str(N(i)),'_p',num2str(p(i)),'.mat']);
load([root_path,data_folder,'data_',num2str(i),'_axis_limit.mat']);
load([root_path,backup_folder,'Q_size',num2str(N(i)),'.mat']);


%% AGM

fprintf('AGM : Press any key to continue.\n');
pause();
clc;

% [~,ng1,x1] = AGM_beta1(d,Q,a,L,x0,lambda,delta,tol,print_output);
% [~,ng2,x2] = AGM_beta2(d,Q,a,L,x0,lambda,delta,tol,mu,print_output);

%% Plotting

% [x,xstar,label] = check_cluster(x1,a,epsilon,handle_outlier);
% 
% title = ['Result : Size = ', num2str(N(i)), ' p = ', num2str(p(i))];
% 
% img_name = ['size',num2str(N(i)),'p',num2str(p(i)),'_result_AGM.png'];
% path = [root_path,img_name];
% 
% [f1] = generate_plot(x,label,path,xlimits,ylimits,title,showfig);

% img_name = ['size',num2str(N(i)),'p',num2str(p(i)),'_result_AGM_xstar.png'];
% path = [root_path,image_folder,img_name];
% [f2] = generate_plot(xstar,label,path,xlimits,ylimits,'test',showfig);

% ax1 = f1.Parent;
% ax2 = f2.Parent;
% fig1 = ax1.Parent;
% fig2 = ax2.Parent;
% spax1 = subplot(1,2,1);
% f1.Parent = spax1;
% spax2 = subplot(1,2,2);
% f2.Parent = spax2;

%% Newton_CG

fprintf('Newton CG : Press any key to continue.\n');
pause();
clc;

options.gamma = 10^(-2);
options.s = 1;
options.sigma = 0.5;
options.maxit = 10;
% tic;
[~,ng,x] = newton_cg(d,Q,a,x0,lambda,delta,tol,options,print_output);
% toc;

%% Plotting

[x,xstar,label] = check_cluster(x,a,epsilon,handle_outlier);

path = [root_path,'size',num2str(N(i)),'p',num2str(p(i)),'_result_NewtonCG.png'];
[~] = generate_plot(x,label,path,xlimits,ylimits,'test',showfig);

% path = [root_path,image_folder,'size',num2str(N(i)),'p',num2str(p(i)),'_result_NewtonCG_xstar.png'];
% [~] = generate_plot(xstar,label,path,xlimits,ylimits,title,showfig);

%% Standard Gradient with Armijo

fprintf('Gradient with Armijo : Press any key to continue.\n');
pause();
clc;

% [~,ng,x] = standard_gradient_armijo(d,Q,a,L,x0,lambda,delta,tol,options,print_output);

%% Plotting

%% Weighted Models

fprintf('Weighted Models ( AGM(beta1) ) : Press any key to continue.\n');
pause();
clc;

% k = [5,10];
theta = 0.5;
k = 5;

[w] = generate_weight(a,k,theta);
L = 1 + N(i)*max(w)/delta;

%% Weighted AGM (beta1)

% [~,ng,x] = AGM_weighted(d,Q,a,L,x0,lambda,delta,tol,w,print_output);

%% Plotting

% [x,xstar,label] = check_cluster(x,a,epsilon,handle_outlier);
% 
% title = ['Result : Size = ', num2str(N(i)), ' p = ', num2str(p(i))];
% img_name = ['size',num2str(N(i)),'p',num2str(p(i)),'_result_AGM_weighted.png'];
% path = [root_path,image_folder,img_name];
% [~] = generate_plot(x,label,path,xlimits,ylimits,title,showfig);
