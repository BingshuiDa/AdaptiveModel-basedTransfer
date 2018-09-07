%This is additional experiment for the CYB reviews. This is to update the
%test problems for multi-objective optimization problem, and to add
%comparison with antuencoder as well. The test problems are WFG problems.

% In WFG test problems: 
% k: the number of position-related parameters
% l: the number of distance-realted parameters

%% Model preparation
% clear 
clc
addpath('moead\')

global maxdim k l M;
f_functions = {'WFG1', 'WFG2', 'WFG3', 'WFG4', 'WFG5', 'WFG6', 'WFG7', 'WFG8', 'WFG9'};
result_TS = {};
alpha = {};
result_NSGA = {};
result_CBR = {};
result_AE = {};
result_MOEAD = {};
M = 2; % number of objectives
k = 4;
l = 20;
maxdim = 24;
dim = k + l;
TrInt = 10;
reps = 30;

for prob = 1:length(f_functions)
    f_source = f_functions;
    f_source(prob) = [];
    f_target = f_functions{prob};

    allmodels_TS = {}; % transfer stacking models
    save('allmodels_TS','allmodels_TS')

    allmodels_CBR = []; % Case-Based reasoning
    save('allmodels_CBR','allmodels_CBR')

    allmodels_AE = {}; % Autoencoder denoising
    save('allmodels_AE', 'allmodels_AE')

    f = testmop(f_source{1},dim);
    L = f.domain(:,1)'; U =  f.domain(:,2)';
    pop = 100;
    no_of_objs = M;
    NSGA_II(f, no_of_objs, L, U, pop);

    f = testmop(f_source{2},dim);
    L = f.domain(:,1)'; U =  f.domain(:,2)';
    pop = 100;
    no_of_objs = M;
    NSGA_II(f, no_of_objs, L, U, pop);

    f = testmop(f_source{3},dim);
    L = f.domain(:,1)'; U =  f.domain(:,2)';
    pop = 100;
    no_of_objs = M;
    NSGA_II(f, no_of_objs, L, U, pop);

    f = testmop(f_source{4},dim);
    L = f.domain(:,1)'; U =  f.domain(:,2)';
    pop = 100;
    no_of_objs = M;
    NSGA_II(f, no_of_objs, L, U, pop);

    f = testmop(f_source{5},dim);
    L = f.domain(:,1)'; U =  f.domain(:,2)';
    pop = 100;
    no_of_objs = M;
    NSGA_II(f, no_of_objs, L, U, pop);

    f = testmop(f_source{6},dim);
    L = f.domain(:,1)'; U =  f.domain(:,2)';
    pop = 100;
    no_of_objs = M;
    NSGA_II(f, no_of_objs, L, U, pop);
    
    f = testmop(f_source{7},dim);
    L = f.domain(:,1)'; U =  f.domain(:,2)';
    pop = 100;
    no_of_objs = M;
    NSGA_II(f, no_of_objs, L, U, pop);
    
    f = testmop(f_source{8},dim);
    L = f.domain(:,1)'; U =  f.domain(:,2)';
    pop = 100;
    no_of_objs = M;
    NSGA_II(f, no_of_objs, L, U, pop);
        
    % Transfer Stacking
    close all

    f = testmop(f_target,dim);
    L = f.domain(:,1)'; U =  f.domain(:,2)';
    pop = 100;
    no_of_objs = M;

    result_TS_prob = [];
    alpha_prob = [];
    result_NSGA_prob = [];
    result_AE_prob = [];
    result_CBR_prob = [];
    result_MOEAD_prob = [];
    for rep = 1:reps
        disp(['Repetition: ', rep])
        disp('Transfer stacking')
        [score_TS, alpha_TS] = TS_NSGA(f, no_of_objs, L, U, pop, true, TrInt);
        result_TS_prob = [result_TS_prob; score_TS];
        alpha_prob = [alpha_prob; alpha_TS];

        disp('No transfer')
        score_NSGA = TS_NSGA(f, no_of_objs, L, U, pop, false);
        result_NSGA_prob = [result_NSGA_prob; score_NSGA];
        
        disp('Autoencoder')
        score_AE = autoencoderEA(f, no_of_objs, L, U, pop, TrInt);
        result_AE_prob = [result_AE_prob; score_AE];
        
        disp('Case-injected')
        score_CBR = CBR_NSGA(f, no_of_objs, L, U, pop);
        result_CBR_prob = [result_CBR_prob; score_CBR];
        
        disp('MOEA/D')
        [pareto1, store1] = moead( f, 'popsize', 100, 'niche', 20, 'iteration', 99, 'method', 'te');
        result_MOEAD_prob = [result_MOEAD_prob; store1];
    end

    result_TS{prob} = result_TS_prob;
    alpha{prob} = alpha_prob;
    result_NSGA{prob} = result_NSGA_prob;
    result_AE{prob} = result_AE_prob;
    result_CBR{prob} = result_CBR_prob;
    result_MOEAD{prob} = result_MOEAD_prob;
end

save('result_WFG.mat', 'result_TS', 'alpha', 'result_NSGA', 'result_CBR', 'result_AE', 'result_MOEAD');
%% 
load result_WFG.mat
f_functions = {'WFG1', 'WFG2', 'WFG3', 'WFG4', 'WFG5', 'WFG6', 'WFG7', 'WFG8', 'WFG9'};
for i = 1:length(f_functions)
    m_TS(i) = mean(result_TS{i}(:,end)); 
    m_NSGA(i) = mean(result_NSGA{i}(:,end)); 
    m_AE(i) = mean(result_AE{i}(:,end)); 
    m_CBR(i) = mean(result_CBR{i}(:,end)); 
    m_MOEAD(i) = mean(result_MOEAD{i}(:,end)); 

    std_TS(i) = std(result_TS{i}(:,end)); 
    std_NSGA(i) = std(result_NSGA{i}(:,end)); 
    std_AE(i) = std(result_AE{i}(:,end)); 
    std_CBR(i) = std(result_CBR{i}(:,end)); 
    std_MOEAD(i) = std(result_MOEAD{i}(:,end)); 
end

%% friedman test
% reps = 30;
% f_functions = {'WFG1', 'WFG2', 'WFG3', 'WFG4', 'WFG5', 'WFG6', 'WFG7', 'WFG8', 'WFG9'};
% for i = 1:length(f_functions)
%     disp(i)
%     m_TS = result_TS{i}(:,end); 
%     m_NSGA = (result_NSGA{i}(:,end)); 
%     m_CBR = (result_CBR{i}(:,end)); 
%     m_MOEAD = (result_MOEAD{i}(:,end)); 
% 
%     p = signrank(m_TS,m_NSGA,'alpha', 0.05);
%     disp(mean(m_TS) - mean(m_NSGA))
%     p = signrank(m_TS,m_CBR,'alpha', 0.05);
%     disp(mean(m_TS) - mean(m_CBR))
%     p = signrank(m_TS,m_MOEAD,'alpha', 0.05);
%     disp(mean(m_TS) - mean(m_MOEAD))
% end


%%
f_functions = {'WFG1', 'WFG2', 'WFG3', 'WFG4', 'WFG5', 'WFG6', 'WFG7'};
for i = 1:length(f_functions)
    m_TS = mean(result_TS{i}); 
    m_NSGA = mean(result_NSGA{i}); 
    m_AE = mean(result_AE{i});
    m_CBR = mean(result_CBR{i}); 
    m_MOEAD = mean(result_MOEAD{i}); 
    x = 100:100:10000;

    figure1 = figure();
    % Create axes
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    plot1 = plot(x(1), m_TS(1), x(1), m_CBR(1), x(1), m_AE(1), x(1), m_NSGA(1), x(1), m_MOEAD(1), ...
    x(1:5:100), m_TS(1:5:100), x(1:5:100), m_CBR(1:5:100), x(1:5:100), m_AE(1:5:100),...
    x(1:5:100), m_NSGA(1:5:100), x(1:5:100), m_MOEAD(1:5:100), ...
    x, m_TS, x, m_CBR, x, m_AE, x, m_NSGA, x, m_MOEAD, ...
        'LineWidth',1,'Parent',axes1);
    set(plot1(1), 'Marker','o', 'MarkerFaceColor',[0 0 0],'Color',[0 0 0], 'LineWidth', 2);
    set(plot1(2), 'Marker','v', 'MarkerFaceColor',[1 0 0],'Color',[1 0 0], 'LineWidth', 2);
    set(plot1(3), 'Marker','+','Color',[0 0 1], 'LineWidth', 2, 'MarkerSize', 8);
    set(plot1(4), 'Marker','square', 'MarkerFaceColor',[0 0 0],'Color',[0.5 0.5 0], 'LineWidth', 2);
    set(plot1(5), 'Marker', '*', 'Color', [0 0.5 0.5], 'LineWidth', 2)
    
    set(plot1(6), 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
    set(plot1(7), 'LineStyle', 'none', 'Marker','v', 'MarkerFaceColor',[1 0 0],'Color',[1 0 0]);
    set(plot1(8), 'LineStyle', 'none', 'Marker','+', 'MarkerSize', 8, 'Color',[0 0 1], 'LineWidth', 2);
    set(plot1(9), 'LineStyle', 'none', 'Marker','square', 'MarkerFaceColor', [0.5, 0.5, 0],'Color',[0.5 0.5 0]);
    set(plot1(10), 'LineStyle', 'none', 'Marker', '*', 'Color', [0 0.5 0.5])
    
    set(plot1(11), 'LineWidth', 2,'Color',[0 0 0]);
    set(plot1(12), 'LineWidth', 2,'Color',[1 0 0]);
    set(plot1(13), 'LineWidth', 2,'Color',[0 0 1]);
    set(plot1(14), 'LineWidth', 2,'Color',[0.5 0.5 0]);
    set(plot1(15), 'LineWidth', 2,'Color',[0 0.5 0.5]);
    set(gca, 'YScale', 'log')

    box(axes1,'on');
    % Create legend
    legend1 = legend('AMTEA', 'TCIEA', 'AEEA', 'NSGA-II', 'MOEA/D');
    set(legend1, 'FontSize', 15, 'Location', 'northeast');
    xlabel('Number of Function Evaluations', 'FontSize', 15)
    ylabel('Averaged IGD Value', 'FontSize', 15)
end
%%
% f_functions = {'ZDT1', 'ZDT2', 'ZDT3', 'ZDT4', 'ZDT6'};
% m_TS = zeros(1, 100);
% m_NSGA = zeros(1, 100);
% m_CBR = zeros(1, 100);
% m_MOEAD = zeros(1, 100);
% for i = 1:length(f_functions)
%     m_TS = m_TS + mean(result_TS{i}); 
%     m_NSGA = m_NSGA+ mean(result_NSGA{i}); 
%     m_CBR = m_CBR + mean(result_CBR{i}); 
%     m_MOEAD = m_MOEAD + mean(result_MOEAD{i}); 
% end
% 
% m_TS = m_TS/5;
% m_NSGA = m_NSGA/5;
% m_CBR = m_CBR/5;
% m_MOEAD = m_MOEAD/5;
% 
% x = 100:100:10000;
% 
% figure1 = figure();
% % Create axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% plot1 = plot(x(1), m_TS(1), x(1), m_CBR(1), x(1), m_NSGA(1), x(1), m_MOEAD(1), ...
%     x(1:5:100), m_TS(1:5:100), x(1:5:100), m_CBR(1:5:100), ...
%     x(1:5:100), m_NSGA(1:5:100), x(1:5:100), m_MOEAD(1:5:100), ...
%     x, m_TS, x, m_CBR, x, m_NSGA, x, m_MOEAD, ...
%         'LineWidth',1,'Parent',axes1);
% set(plot1(1), 'Marker','o', 'MarkerFaceColor',[0 0 0],'Color',[0 0 0], 'LineWidth', 2);
% set(plot1(2), 'Marker','v', 'MarkerFaceColor',[1 0 0],'Color',[1 0 0], 'LineWidth', 2);
% set(plot1(3), 'Marker','+','Color',[0 0 1], 'LineWidth', 2);
% set(plot1(4), 'Marker','square', 'MarkerFaceColor',[0 0 0],'Color',[0.5 0.5 0], 'LineWidth', 2);
% set(plot1(5), 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
% set(plot1(6), 'LineStyle', 'none', 'Marker','v', 'MarkerFaceColor',[1 0 0],'Color',[1 0 0]);
% set(plot1(7), 'LineStyle', 'none', 'Marker','+', 'MarkerSize', 8, 'Color',[0 0 1], 'LineWidth', 2);
% set(plot1(8), 'LineStyle', 'none', 'Marker','square', 'MarkerFaceColor', [0.5, 0.5, 0],'Color',[0.5 0.5 0]);
% set(plot1(9), 'LineWidth', 2,'Color',[0 0 0]);
% set(plot1(10), 'LineWidth', 2,'Color',[1 0 0]);
% set(plot1(11), 'LineWidth', 2,'Color',[0 0 1]);
% set(plot1(12), 'LineWidth', 2,'Color',[0.5 0.5 0]);
% set(gca, 'YScale', 'log')
% 
% box(axes1,'on');
% % Create legend
% legend1 = legend('AMTEA', 'TCIEA', 'NSGA-II', 'MOEA/D');
% set(legend1, 'FontSize', 15, 'Location', 'northeast');
% xlabel('Number of Function Evaluations', 'FontSize', 15)
% ylabel('Averaged IGD Value', 'FontSize', 15)

%%
reps = 30;
ind = [1,2,3,4,5,6,7,8,9];
num_alpha = 10;
n_tasks = length(f_functions);
for i = 1:length(ind)
    alpha_prob = alpha{ind(i)};
    alpha_Matrix = zeros(num_alpha,n_tasks);
    for rep = 1:reps
        alpha_Matrix = alpha_Matrix + alpha_prob((rep-1)*num_alpha+1:rep*num_alpha,:);
    end
    alpha_Matrix = alpha_Matrix/reps;

    error_Matrix = zeros(num_alpha, n_tasks);
    for rep = 1:rep
    	error_Matrix = error_Matrix + (alpha_prob((rep-1)*num_alpha+1:rep*num_alpha,:) - alpha_Matrix).^2;
    end
    error_Matrix = sqrt(error_Matrix/reps);

    figure2 = figure;

    % Create axes
    axes2 = axes('Parent',figure2);
    hold(axes2,'on');
    x = 2:num_alpha:100;
    plot2 = plot(x, alpha_Matrix(:,1),x, alpha_Matrix(:,2),...
        x, alpha_Matrix(:,3),x, alpha_Matrix(:,4), ...
        x, alpha_Matrix(:,5), x, alpha_Matrix(:,6),'Parent',axes2);
    set(plot2(1),'LineWidth',2,'Color',[0 0 1],'Marker','o', 'MarkerFaceColor',[0 0 1]);
    set(plot2(2),'LineWidth',2,'Color',[1 0 0],'Marker','v', 'MarkerFaceColor', [1 0 0]);
    set(plot2(3),'LineWidth',2,'Color',[0.5 0.5 0],'Marker','+', 'MarkerSize', 8);
    set(plot2(4),'LineWidth',2,'Color',[0 0 0],'Marker','square', 'MarkerFaceColor', [0 0 0]);
    set(plot2(5),'LineWidth',2,'Color',[0 0.5 0.5],'Marker','*', 'MarkerSize', 8);
    set(plot2(6),'LineWidth',2,'Color',[0.3 0.3 0.3],'Marker','x', 'MarkerFaceColor', [0 0 0]);
    
    shadedErrorBar(x, alpha_Matrix(:, 1)', error_Matrix(:, 1)','lineprops',{'Color',[0 0 1]},'transparent',1);
    shadedErrorBar(x, alpha_Matrix(:, 2)', error_Matrix(:, 2)','lineprops',{'Color',[1 0 0]},'transparent',1);
    shadedErrorBar(x, alpha_Matrix(:, 3)', error_Matrix(:, 3)','lineprops',{'Color',[0.5 0.5 0]},'transparent',1);
    shadedErrorBar(x, alpha_Matrix(:, 4)', error_Matrix(:, 4)','lineprops',{'Color',[0 0 0]},'transparent',1);
    shadedErrorBar(x, alpha_Matrix(:, 5)', error_Matrix(:, 5)','lineprops',{'Color',[0 0.5 0.5]},'transparent',1);
    shadedErrorBar(x, alpha_Matrix(:, 6)', error_Matrix(:, 6)','lineprops',{'Color',[0.3 0.3 0.3]},'transparent',1);

    ylim([0,Inf])
    if i == 1
        legend2 = legend('WFG2 (source 1)', 'WFG3 (source 2)', 'WFG4 (source 3)', 'WFG5 (source 4)', 'WFG6 (source 5)', 'WFG7 (source 6)');
    elseif i == 2
        legend2 = legend('WFG1 (source 1)', 'WFG3 (source 2)', 'WFG4 (source 3)', 'WFG5 (source 4)', 'WFG6 (source 5)', 'WFG7 (source 6)');
    elseif i == 3
        legend2 = legend('WFG1 (source 1)', 'WFG2 (source 2)', 'WFG4 (source 3)', 'WFG5 (source 4)', 'WFG6 (source 5)', 'WFG7 (source 6)');
    elseif i == 4
        legend2 = legend('WFG1 (source 1)', 'WFG2 (source 2)', 'WFG3 (source 3)', 'WFG5 (source 4)', 'WFG6 (source 5)', 'WFG7 (source 6)');
    elseif i == 5
        legend2 = legend('WFG1 (source 1)', 'WFG2 (source 2)', 'WFG3 (source 3)', 'WFG4 (source 4)', 'WFG6 (source 5)', 'WFG7 (source 6)');
    elseif i == 6
        legend2 = legend('WFG1 (source 1)', 'WFG2 (source 2)', 'WFG3 (source 3)', 'WFG4 (source 4)', 'WFG5 (source 5)', 'WFG7 (source 6)');
    else
        legend2 = legend('WFG1 (source 1)', 'WFG2 (source 2)', 'WFG3 (source 3)', 'WFG4 (source 4)', 'WFG5 (source 5)', 'WFG6 (source 6)');
    end
    set(legend2, 'FontSize', 15);
    xlabel('Number of Generations', 'FontSize', 15)
    ylabel('Transfer Coefficients (\alpha''s)', 'FontSize', 15)
end