clc;
clear;
close all;
format long;

%%

root_path = './';           % Root path
image_folder = 'plots/';    % Image folder
data_folder = 'datasets/';  % Dataset folder
backup_folder = 'backup/';  % Backup folder for auxiliary matrices
result_folder = 'results/'; % Results folder

% print_output = false;       % Print output
print_output = true;
% showfig = true;             % Figure display
showfig = false;

%%

handle_outlier = false;     % Outlier checking and removing

lambdas = 0.05:0.005:0.1;   % 11 lambdas

% tol = 10^(-2);              % Stopping criterion
tol = 10^(-1);              % Stopping criterion
delta = 10^(-4);            % Huber norm

options.gamma = 0.1;        % Backtracking
options.s = 1;              % Backtracking
options.sigma = 0.5;        % Backtracking

options.maxit = 10;         % Maximum number of CG iterations

kn = 5;                      % knn
theta = 0.5;                % knn

epsilon = 10^(-2);          % Check cluster

%% Datasets

p = [2,3,4];
n = [200,600,400];
Num = size(p,2);            % Num of datasets

%% 

skip = [1,2];
SkipMethod = [4];

for i = 1:Num
    if ismember(i,skip)
        continue;
    end
    dataset_name = ['data_',num2str(i),'_size',num2str(n(i)),'_p',num2str(p(i)),'.mat'];
    auxmatrix_name = ['Q_size',num2str(n(i)),'.mat'];
    load([root_path,data_folder,dataset_name]);
    load([root_path,backup_folder,auxmatrix_name]);
    
    % Initialization
    N = size(a,1);
    d = size(a,2);
    x0 = sparse(N,d);
    mu = 1;
    
    for j = 1:size(lambdas,2)
        lambda = lambdas(j);
        L = 1 + N*lambda/delta;
        
        k = 1;
        if ~ismember(k,SkipMethod)
            % AGM beta1
            tic;
            [iter1,ng1,x1] = AGM_beta1(d,Q,a,L,x0,lambda,delta,tol,print_output);
            t1 = toc;
        end
        k = k + 1;
        if ~ismember(k,SkipMethod)
            % AGM beta2
            tic;
            [iter2,ng2,x2] = AGM_beta2(d,Q,a,L,x0,lambda,delta,tol,mu,print_output);
            t2 = toc;
        end
        k = k + 1;
        if ~ismember(k,SkipMethod)
            % NewtonCG
            tic;
            [iter3,ng3,x3] = newton_cg(d,Q,a,x0,lambda,delta,tol,options,print_output);
            t3 = toc;
        end
        k = k + 1;
        if ~ismember(k,SkipMethod)
            % Gradient with Armijo
            tic;
            [iter4,ng4,x4] = standard_gradient_armijo(d,Q,a,L,x0,lambda,delta,tol,options,print_output);
            t4 = toc;
        end
        k = k + 1;
        if ~ismember(k,SkipMethod)
            % AGM beta1 weighted
            tic;
            [w] = generate_weight(a,kn,theta);
            L = 1 + N*max(w)/delta;
            [iter5,ng5,x5] = AGM_weighted(d,Q,a,L,x0,lambda,delta,tol,w,print_output);
            t5 = toc;
        end
        
        % Save results to .mat files
        filename = ['result_data_',num2str(i),'_lambda',num2str(j),'_',num2str(lambda),'.mat'];
%         save([root_path,result_folder,filename],'iter1','ng1','x1','t1','iter2','ng2','x2','t2',...
%             'iter3','ng3','x3','t3','iter4','ng4','x4','t4','iter5','ng5','x5','t5');
%         save([root_path,filename],'iter3','ng3','x3','t3');
%         save([root_path,result_folder,filename],'iter1','ng1','x1','t1','iter2','ng2','x2','t2',...
%             'iter3','ng3','x3','t3','iter5','ng5','x5','t5');
    end
end

%% Plotting

skip = [1,2];

plot_xstar = false;

for i = 1:Num
    if ismember(i,skip)
        continue;
    end
    for j = 1:size(lambdas,2)
        lambda = lambdas(j);
        
        dataset_name = ['data_',num2str(i),'_size',num2str(n(i)),'_p',num2str(p(i)),'.mat'];
        filename = ['result_data_',num2str(i),'_lambda',num2str(j),'_',num2str(lambda),'.mat'];
        misc = ['data_',num2str(i),'_axis_limit.mat'];
        load([root_path,data_folder,dataset_name]);
        load([root_path,result_folder,filename]);
        load([root_path,data_folder,misc]);
        
        image_title = ['Size = ',num2str(n(i)),' ','p = ',num2str(p(i)),' ','lambda = ',num2str(lambda),' ','ElapsedTime = '];
        img_name = ['data_',num2str(i),'_results_'];
        path = [root_path,image_folder,img_name];
        
        % AGM beta1
        [x1,xstar1,label1] = check_cluster(x1,a,epsilon,handle_outlier);
        [f1] = generate_plot(x1,label1,[path,'AGMbeta1_lambda',num2str(j),'_',num2str(lambda),'.png'],xlimits,ylimits,{'AGMbeta1',[image_title,num2str(t1)]},showfig);
        % AGM beta2
        [x2,xstar2,label2] = check_cluster(x2,a,epsilon,handle_outlier);
        [f2] = generate_plot(x2,label2,[path,'AGMbeta2_lambda',num2str(j),'_',num2str(lambda),'.png'],xlimits,ylimits,{'AGMbeta2',[image_title,num2str(t2)]},showfig);
        % NewtonCG
        [x3,xstar3,label3] = check_cluster(x3,a,epsilon,handle_outlier);
        [f3] = generate_plot(x3,label3,[path,'NewtonCG_lambda',num2str(j),'_',num2str(lambda),'.png'],xlimits,ylimits,{'NewtonCG',[image_title,num2str(t3)]},showfig);
%         % Gradient with Armijo
%         [x4,xstar4,label4] = check_cluster(x4,a,epsilon,handle_outlier);
%         [f4] = generate_plot(x4,label4,[path,'GM_lambda',num2str(j),'_',num2str(lambda),'.png'],xlimits,ylimits,{'GM',[image_title,num2str(t4)]},showfig);
        % AGM beta1 weighted
        [x5,xstar5,label5] = check_cluster(x5,a,epsilon,handle_outlier);
        [f5] = generate_plot(x5,label5,[path,'AGMbeta1_weighted_lambda',num2str(j),'_',num2str(lambda),'.png'],xlimits,ylimits,{'AGMbeta1Weighted',[image_title,num2str(t5)]},showfig);
        
        if plot_xstar == true
            [f11] = generate_plot(xstar1,label,[path,'AGMbeta1_xstar_lambda',num2str(j),'_',num2str(lambda),'.png'],xlimits,ylimits,{'AGMbeta1(xstar)',[image_title,num2str(t1)]},showfig);
            [f22] = generate_plot(xstar2,label,[path,'AGMbeta2_xstar_lambda',num2str(j),'_',num2str(lambda),'.png'],xlimits,ylimits,{'AGMbeta2(xstar)',[image_title,num2str(t2)]},showfig);
            [f33] = generate_plot(xstar3,label,[path,'NewtonCG_xstar_lambda',num2str(j),'_',num2str(lambda),'.png'],xlimits,ylimits,{'NewtonCG(xstar)',[image_title,num2str(t3)]},showfig);
%             [f44] = generate_plot(xstar4,label,[path,'GM_xstar_lambda',num2str(j),'_',num2str(lambda),'.png'],xlimits,ylimits,{'GM(xstar)',[image_title,num2str(t4)]},showfig);
            [f55] = generate_plot(xstar5,label,[path,'AGMbeta1_weighted_xstar_lambda',num2str(j),'_',num2str(lambda),'.png'],xlimits,ylimits,{'AGMbeta1Weighted(xstar)',[image_title,num2str(t5)]},showfig);
        end
        
        % Convergence Plot
        f = figure;
        if showfig == false
            set(f,'visible','off');
        end
        hold on;
        plot(1:size(ng1,1),log10(ng1),'LineWidth',1);
        plot(1:size(ng2,1),log10(ng2),'LineWidth',1);
        plot(1:size(ng3,1),log10(ng3),'LineWidth',1);
%         plot(1:size(ng4,1),log10(ng4),'LineWidth',1);
        plot(1:size(ng5,1),log10(ng5),'LineWidth',1);
        hold off;
        ylim([-1 inf]);
        xlim([0 inf]);
        xlabel('Number of iterations');
        title({'Convergence Plot',['lambda = ',num2str(lambda)]});
%         legend({'AGMbeta1','AGMbeta2','NewtonCG','GM','AGMbeta1-weighted'},'Location','southeast');
        legend({'AGMbeta1','AGMbeta2','NewtonCG','AGMbeta1-weighted'},'Location','southeast');
        saveas(gcf,[root_path,image_folder,'data_',num2str(i),'_convergence_plot_lambda',num2str(j),'_',num2str(lambda),'.png']);
    end
end

%% Extensions

clc;
close all;

load('./datasets/extension/grad_norm_1_adam.mat');
g1 = a';
load('./datasets/extension/grad_norm_1_agm_l1.mat');
g2 = a';
load('./datasets/extension/grad_norm_1_agm_linf.mat');
g3 = a';
load('./datasets/extension/grad_norm_1_agm_unconvex.mat');
g4 = a';
load('./datasets/extension/grad_norm_1_sgd.mat');
g5 = a';

load([root_path,result_folder,'result_data_1_lambda11_0.1.mat']);
f = figure;
if showfig == false
    set(f,'visible','off');
end
hold on;
plot(1:size(ng1,1),log10(ng1),'LineWidth',1);
plot(1:size(g2,1),log10(g2),'LineWidth',1);
plot(1:size(g3,1),log10(g3),'LineWidth',1);
plot(1:size(g4,1),log10(g4),'LineWidth',1);
hold off;
ylim([-1 inf]);
xlim([1 10000]);
xlabel('Number of iterations');
title({'Convergence Plot','lambda = 0.1'});
legend({'AGMbeta1','AGM-L1','AGM-Linf','AGM-Unconvex'},'Location','southeast');
% saveas(gcf,[root_path,image_folder,'data_1_convergence_plot_lambda_0.1_extension1.png']);


f = figure;
if showfig == false
    set(f,'visible','off');
end
hold on;
plot(1:size(ng1,1),log10(ng1),'LineWidth',1);
plot(1:size(ng3,1),log10(ng3),'LineWidth',1);
plot(1:size(g1,1),log10(g1),'LineWidth',1);
plot(1:size(g5,1),log10(g5),'LineWidth',1);
hold off;
ylim([-1 inf]);
xlim([1 10000]);
xlabel('Number of iterations');
title({'Convergence Plot','lambda = 0.1'});
legend({'AGMbeta1','NewtonCG','ADAM','SGD'},'Location','southeast');
% saveas(gcf,[root_path,image_folder,'data_1_convergence_plot_lambda_0.1_extension2.png']);
