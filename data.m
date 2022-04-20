clc;
clear;
close all;
format long;

%%

root_path = './';           % Root path
image_folder = 'plots/';    % Image folder
data_folder = 'datasets/';  % Dataset folder
backup_folder = 'backup/';  % Backup folder for auxiliary matrices

margin = 5;                 % Margin for plotting

Num = 4;                    % Num of datasets
d   = 2;                    % Dimension of data points
p   = [2,3,4,10];           % Num of clusters
N   = [200,600,400,3000];  % Num of observations
n   = N./p;                 % Num of data points for each cluster

c1 = [0,0;15,15];
c2 = [0,0;15,30;30,5];
c3 = [0,0;10,20;20,0;30,20];
c4 = [15,10;10,5;5,10;15,0;0,5;20,5;0,15;20,15;10,20;10,15];
c = {c1,c2,c3,c4};

sigma1 = [2.5;2.5];
sigma2 = [5;5;5];
sigma3 = [5;5;5;5];
sigma4 = [1.7;2;1.3;2;1.5;2;2;2.1;1.7;1.5];
sigma = {sigma1,sigma2,sigma3,sigma4};


%% Data Generation

showfig = true;
% showfig = false;

Skip = [1,2,3,4];

for i = 1:Num
    if ismember(i,Skip)
        continue;
    end
    centroids = cell2mat(c(i));
    sigmas = cell2mat(sigma(i));
    [a] = generate_data(d,p(i),n(i),centroids,sigmas);
    save([root_path,data_folder,'data_',num2str(i),'_size',num2str(N(i)),'_p',num2str(p(i))],'a');
    fprintf('%d : Data Generated :\tsize : %d\t Num of clusters : %d\n',i,N(i),p(i));
    
    f = figure;  
    if showfig == false
        set(f,'visible','off');
    end
    hold on;
    for j = 1:p(i)
        start = (j-1)*n(i)+1;
        stop = j*n(i);
        scatter(a(start:stop,1),a(start:stop,2),[],j.*ones(1,stop-start+1),'filled');
    end
    title(['Size = ', num2str(N(i)), ' p = ', num2str(p(i))]);
    [xmin,indmin] = min(centroids(:,1));
    [xmax,indmax] = max(centroids(:,1));
    xlimits = [xmin - 3*sigmas(indmin) - margin, xmax + 3*sigmas(indmax) + margin];
    xlim(xlimits);
    [ymin,indmin] = min(centroids(:,2));
    [ymax,indmax] = max(centroids(:,2));
    ylimits = [ymin - 3*sigmas(indmin) - margin, ymax + 3*sigmas(indmax) + margin];
    ylim(ylimits);
    hold off;
    saveas(gcf,[root_path,image_folder,'data_',num2str(i),'_size',num2str(N(i)),'_p',num2str(p(i)),'.png']);
    save([root_path,data_folder,'data_',num2str(i),'_axis_limit.mat'],'xlimits','ylimits');
end

%% Auxiliary Matrices Generation

Skip = [1,2,3,4];

for i = 1:Num
    if ismember(i,Skip)
        continue;
    end
    [Q] = generate_Q(N(i));
    save([root_path,backup_folder,'Q_size',num2str(N(i)),'.mat'],'Q');
end