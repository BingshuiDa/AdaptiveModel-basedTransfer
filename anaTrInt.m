%This is additional experiment for the CYB reviews. This is to study the
%effect of transfer interval. Compare the different experimental result
%when TrInt=2,5,10

%%
clear
clc
allmodels_TS = {};
save('allmodels_TS','allmodels_TS');
allmodels_CBR = [];
save('allmodels_CBR', 'allmodels_CBR');

reps = 30;
TrInt = [2, 5, 10];

BGA('onemax',100,100);
BGA('onemin',100,100);
result_TS = {};
alpha = {};
result_CBR = {};
trans.transfer = false;
result_GA = TSBGA('trap5',100,reps,trans);

for i = 1:length(TrInt)
    trans.transfer = true;
    trans.TrInt = TrInt(i);
    [result_TS_prob, alpha_prob] = TSBGA('trap5',100,reps,trans);
    result_TS = [result_TS; result_TS_prob];
    alpha = [alpha; alpha_prob];
    
    trans.transfer = false;
    result_CBR_prob = CBRBGA('trap5',100,reps,TrInt(i));
    result_CBR = [result_CBR; result_CBR_prob];
end
save('result_trap_TrInt.mat', 'result_TS', 'alpha', 'result_GA', 'result_CBR');

%%
% clc
% allmodels_TS = {};
% save('allmodels_TS','allmodels_TS');
% allmodels_CBR = [];
% save('allmodels_CBR', 'allmodels_CBR');
% 
% reps = 30;
% TrInt = [2, 5, 10];
% 
% BGA('onemax',100,100);
% trans.transfer = true;
% trans.transfer = 2;
% [~, ~] = TSBGA('trap5',100, 1, trans, true, 100);
% BGA('trap5',100,80);
% 
% result_TS = {};
% alpha = {};
% result_CBR = {};
% trans.transfer = false;
% result_GA = TSBGA('onemin',100,reps,trans);
% 
% for i = 1:length(TrInt)
%     trans.transfer = true;
%     trans.TrInt = TrInt(i);
%     [result_TS_prob, alpha_prob] = TSBGA('onemin',100,reps, trans);
%     result_TS = [result_TS; result_TS_prob];
%     alpha = [alpha; alpha_prob];
%     
%     trans.transfer = false;
%     result_CBR_prob = CBRBGA('onemin',100,reps,TrInt(i));
%     result_CBR = [result_CBR; result_CBR_prob];
% end
% 
% save('result_onemin_TrInt.mat', 'result_TS', 'alpha', 'result_GA', 'result_CBR')


%% draw figure
m_TS_1 = mean(result_TS{1}); 
m_TS_2 = mean(result_TS{2}); 
m_TS_3 = mean(result_TS{3}); 
m_GA = mean(result_GA); 
m_CBR_1 = mean(result_CBR{1}); 
m_CBR_2 = mean(result_CBR{2}); 
m_CBR_3 = mean(result_CBR{3}); 
% m_DT = mean(result_DT);

figure1 = figure;

x = 50:50:5000;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

plot1 = plot(x(1), m_TS_1(1), x(1), m_TS_2(1), x(1), m_TS_3(1), x(1), m_GA(1), x(1), m_CBR_1(1), ...
      x(1), m_CBR_2(1), x(1), m_CBR_3(1), x(1:5:100), m_TS_1(1:5:100), x(1:5:100), m_TS_2(1:5:100), x(1:5:100), m_TS_3(1:5:100), x(1:5:100), m_GA(1:5:100), ...
      x(1:5:100), m_CBR_1(1:5:100), x(1:5:100), m_CBR_2(1:5:100), x(1:5:100), m_CBR_3(1:5:100), ...
      x, m_TS_1, x, m_TS_2, x, m_TS_3, x, m_GA, x, m_CBR_1, x, m_CBR_2, x, m_CBR_3, ...
          'LineWidth',1,'Parent',axes1);
set(plot1(1), 'Marker', 'o', 'MarkerFaceColor',[0 0 0],'Color',[0 0 0], 'LineWidth', 2);
set(plot1(2), 'Marker', '*', 'MarkerSize', 8, 'MarkerFaceColor',[0 0.5 0],'Color',[0 0.5 0], 'LineWidth', 2);
set(plot1(3), 'Marker', '^', 'MarkerFaceColor',[0 0 0.5],'Color',[0 0 0.5], 'LineWidth', 2);
set(plot1(4), 'Marker','v', 'MarkerFaceColor',[1 0 0], 'Color',[1 0 0], 'LineWidth', 2);
set(plot1(5), 'Marker','+', 'MarkerSize', 8, 'Color',[0 0 1], 'LineWidth', 2);
set(plot1(6), 'Marker','x', 'MarkerSize', 8, 'Color',[0 0.5 1], 'LineWidth', 2);
set(plot1(7), 'Marker','square', 'MarkerSize', 8, 'Color',[0.5 0 1], 'LineWidth', 2);

set(plot1(8), 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor',[0 0 0], 'Color',[0 0 0]);
set(plot1(9), 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 8, 'MarkerFaceColor',[0 0.5 0], 'Color',[0 0.5 0]);
set(plot1(10), 'LineStyle', 'none', 'Marker', '^', 'MarkerFaceColor',[0 0 0.5], 'Color',[0 0 0.5]);
set(plot1(11), 'LineStyle', 'none', 'Marker','v', 'MarkerFaceColor',[1 0 0],'Color',[1 0 0]);
set(plot1(12), 'LineStyle', 'none','Color',[0 0 1], 'Marker','+', 'MarkerSize', 8, 'LineWidth', 2);
set(plot1(13), 'LineStyle', 'none','Color',[0 0.5 1], 'Marker','x', 'MarkerSize', 8, 'LineWidth', 2);
set(plot1(14), 'LineStyle', 'none','Color',[0.5 0 1], 'Marker','square', 'MarkerSize', 8, 'LineWidth', 2);

set(plot1(15), 'LineWidth', 2,'Color',[0 0 0]);
set(plot1(16), 'LineWidth', 2,'Color',[0 0.5 0]);
set(plot1(17), 'LineWidth', 2,'Color',[0 0 0.5]);
set(plot1(18), 'LineWidth', 2,'Color',[1 0 0]);
set(plot1(19), 'LineWidth', 2,'Color',[0 0 1]);
set(plot1(20), 'LineWidth', 2,'Color',[0 0.5 1]);
set(plot1(21), 'LineWidth', 2,'Color',[0.5 0 1]);

% ylim([80, 100])
box(axes1,'on');
% Create legend
legend1 = legend('AMTEA $\Delta=2$', 'AMTEA $\Delta=5$', 'AMTEA $\Delta=10$', ...
    'CEA', 'TCIEA $\Delta=2$', 'TCIEA $\Delta=5$', 'TCIEA $\Delta=10$');
set(legend1, 'FontSize', 15, 'Location', 'southeast', 'interpreter','latex');
xlabel('Number of Function Evaluations', 'FontSize', 15, 'interpreter','latex')
ylabel('Averaged Objective Value', 'FontSize', 15, 'interpreter','latex')

%%
load result_trap_TrInt.mat
reps = 30;
alpha_Matrix = zeros(49,3);
for rep = 1:reps
  alpha_Matrix = alpha_Matrix + alpha{1, rep};
end
alpha_Matrix = alpha_Matrix/reps;

error_Matrix = zeros(49,3);
for rep = 1:reps
    error_Matrix = error_Matrix + (alpha{1, rep} - alpha_Matrix).^2;
end
error_Matrix = sqrt(error_Matrix/reps);

figure2 = figure;

% Create axes
axes2 = axes('Parent',figure2);
hold(axes2,'on');
x = 2:2:99;
plot2 = plot(x, alpha_Matrix(:,1), x, alpha_Matrix(:,2), 'Parent',axes2);
set(plot2(1), 'LineWidth', 2, 'Color',[0 0 1],'Marker','o', 'MarkerFaceColor', [0, 0, 1]);
set(plot2(2), 'LineWidth', 2, 'Color',[1 0 0],'Marker','+', 'MarkerSize', 8);

shadedErrorBar(x,alpha_Matrix(:, 1)',0.3*error_Matrix(:, 1)','lineprops',{'Color',[0 0 1]},'transparent',1);
shadedErrorBar(x,alpha_Matrix(:, 2)',0.3*error_Matrix(:, 2)','lineprops',{'Color',[1 0 0]},'transparent',1);

ylim([0,0.6])
box(axes2,'on');
% Create legend
legend2 = legend( 'one-max (source 1)', 'one-min (source 2)');
set(legend2, 'FontSize', 15, 'interpreter','latex');
xlabel('Number of Generations', 'FontSize', 15, 'interpreter','latex')
ylabel('Transfer Coefficients ($\alpha_k$''s)', 'FontSize', 15, 'interpreter','latex')