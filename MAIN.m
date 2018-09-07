% clear 
clc
allmodels_TS = {};
save('allmodels_TS','allmodels_TS');
allmodels_CBR = [];
save('allmodels_CBR', 'allmodels_CBR');

reps = 30;
TrInt = 2;

BGA('onemax',100,100);
BGA('onemin',100,100);
trans.transfer = true;
trans.TrInt = TrInt;
[result_TS, alpha] = TSBGA('trap5',100,reps, trans);
trans.transfer = false;
result_GA = TSBGA('trap5',100,reps,trans);
result_CBR = CBRBGA('trap5',100,reps,TrInt);

%%
clc
allmodels_TS = {};
save('allmodels_TS','allmodels_TS');
allmodels_CBR = [];
save('allmodels_CBR', 'allmodels_CBR');

reps = 30;
TrInt = 2;

BGA('onemax',100,100);
trans.transfer = true;
trans.transfer = TrInt;
[~, ~] = TSBGA('trap5',100, 1, trans, true, 100);
BGA('trap5',100,80);

trans.transfer = true;
trans.transfer = TrInt;
[result_TS, alpha] = TSBGA('onemin',100,reps,trans);
trans.transfer = false;
result_GA = TSBGA('onemin',100,reps,trans);
result_CBR = CBRBGA('onemin',100,reps,TrInt);

save('result_onemin.mat', 'result_TS', 'alpha', 'result_GA', 'result_CBR')

%% draw figure
m_TS = mean(result_TS); 
m_GA = mean(result_GA); 
m_CBR = mean(result_CBR); 
% m_DT = mean(result_DT);

figure1 = figure;

x = 50:50:5000;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

plot1 = plot(x(1), m_TS(1), x(1), m_GA(1), x(1), m_CBR(1), ...
      x(1:5:100), m_TS(1:5:100), x(1:5:100), m_GA(1:5:100), ...
      x(1:5:100), m_CBR(1:5:100), ...
      x, m_TS, x, m_GA, x, m_CBR, ...
          'LineWidth',1,'Parent',axes1);
set(plot1(1), 'Marker', 'o', 'MarkerFaceColor',[0 0 0],'Color',[0 0 0], 'LineWidth', 2);
set(plot1(2), 'Marker','v', 'MarkerFaceColor',[1 0 0], 'Color',[1 0 0], 'LineWidth', 2);
set(plot1(3), 'Marker','+', 'MarkerSize', 8, 'Color',[0 0 1], 'LineWidth', 2);
set(plot1(4), 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor',[0 0 0], 'Color',[0 0 0]);
set(plot1(5), 'LineStyle', 'none', 'Marker','v', 'MarkerFaceColor',[1 0 0],'Color',[1 0 0]);
set(plot1(6), 'LineStyle', 'none','Color',[0 0 1], 'Marker','+', 'MarkerSize', 8, 'LineWidth', 2);
set(plot1(7), 'LineWidth', 2,'Color',[0 0 0]);
set(plot1(8), 'LineWidth', 2,'Color',[1 0 0]);
set(plot1(9), 'LineWidth', 2,'Color',[0 0 1]);

ylim([60, 101])
box(axes1,'on');
% Create legend
legend1 = legend('AMTEA', 'CEA', 'TCIEA');
set(legend1, 'FontSize', 15, 'Location', 'southeast', 'interpreter','latex');
xlabel('Number of Function Evaluations', 'FontSize', 15, 'interpreter','latex')
ylabel('Averaged Objective Value', 'FontSize', 15, 'interpreter','latex')

%%
load result_onemin.mat
reps = 30;
alpha_Matrix = zeros(49,4);
for rep = 1:reps
  alpha_Matrix = alpha_Matrix + alpha{rep};
end
alpha_Matrix = alpha_Matrix/reps;

error_Matrix = zeros(49,4);
for rep = 1:reps
    error_Matrix = error_Matrix + (alpha{rep} - alpha_Matrix).^2;
end
error_Matrix = sqrt(error_Matrix/reps);

figure2 = figure;

% Create axes
axes2 = axes('Parent',figure2);
hold(axes2,'on');
x = 2:2:99;
plot2 = plot(x, alpha_Matrix(:,1), x, alpha_Matrix(:,2), x, alpha_Matrix(:,3), 'Parent',axes2);
set(plot2(1), 'LineWidth', 2, 'Color',[0 0 1],'Marker','o', 'MarkerFaceColor', [0, 0, 1]);
set(plot2(2), 'LineWidth', 2, 'Color',[1 0 0],'Marker','^', 'MarkerFaceColor', [1, 0, 0]);
set(plot2(3), 'LineWidth', 2, 'Color',[0 0 0],'Marker','+','Marker','+', 'MarkerSize', 8);

% errorbar(2:2:100, alpha_Matrix(:, 1)', 0.3*error_Matrix(:, 1)','LineWidth',2,'Color',[0 0 1]);
% errorbar(2:2:100, alpha_Matrix(:, 2)', 0.3*error_Matrix(:, 2)','LineWidth',2,'Color',[1 0 0]);
% errorbar(2:2:99, alpha_Matrix(:, 3)', 0.5*error_Matrix(:, 3)','LineWidth',2,'Color',[0 0 0]);

shadedErrorBar(x,alpha_Matrix(:, 1)',error_Matrix(:, 1)','lineprops',{'Color',[0 0 1]},'transparent',1);
shadedErrorBar(x,alpha_Matrix(:, 2)',error_Matrix(:, 2)','lineprops',{'Color',[1 0 0]},'transparent',1);
shadedErrorBar(x,alpha_Matrix(:, 3)',error_Matrix(:, 3)','lineprops',{'Color',[0 0 0]},'transparent',1);

ylim([0,0.5])
box(axes2,'on');
% Create legend
% legend2 = legend( 'one-max', 'one-min');
legend2 = legend( 'one-max (source 1)', 'trap-5 -- global optimum (source 2)', 'trap-5 -- deceptive local optimum (source 3)');
set(legend2, 'FontSize', 13, 'interpreter','latex');
xlabel('Number of Generations', 'FontSize', 15, 'interpreter','latex')
ylabel('Transfer Coefficients ($\alpha_k$''s)', 'FontSize', 15, 'interpreter','latex')

%%
load result_trap_max.mat
reps = 30;
alpha_Matrix = zeros(49,3);
for rep = 1:reps
  alpha_Matrix = alpha_Matrix + alpha{rep};
end
alpha_Matrix = alpha_Matrix/reps;

error_Matrix = zeros(49,3);
for rep = 1:reps
    error_Matrix = error_Matrix + (alpha{rep} - alpha_Matrix).^2;
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
set(legend2, 'FontSize', 15);
xlabel('Number of Generations', 'FontSize', 15)
ylabel('Transfer Coefficients (\alpha_k''s)', 'FontSize', 15)
