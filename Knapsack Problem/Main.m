clear
clc

load KP_sc_ak
load KP_uc_ak
load KP_wc_ak
load KP_wc_rk
load KP_sc_rk
load KP_uc_rk

allmodels_TS = {};
save('allmodels_TS','allmodels_TS')
allmodels_CBR = [];
save('allmodels_CBR', 'allmodels_CBR')

result_TS = {};
alpha = {};
result_GA = {};
result_CBR = {};

% BGA('onemax',1000,1000);
% BGA('onemin',1000,1000);
KP_BGA(KP_uc_rk,1000,true);
KP_BGA(KP_sc_rk,1000,true);
KP_BGA(KP_wc_rk,1000,true);
KP_BGA(KP_sc_ak,1000,true);

reps = 30;
% KP_wc_ak
[result_TS_prob, alpha_prob] = KP_TSBGA(KP_wc_ak,1000,true,reps);
[result_TS_random_prob, ~] = KP_TSBGA(KP_wc_ak,1000,true,reps,true);% random coefficient
result_GA_prob = KP_TSBGA(KP_wc_ak,1000,false,reps);
result_CBR_prob = KP_CBRGA(KP_wc_ak,1000,reps);
result_TS{1} = result_TS_prob;
alpha{1} = alpha_prob;
result_TS_random{1} = result_TS_random_prob;
result_GA{1} = result_GA_prob;
result_CBR{1} = result_CBR_prob;

% KP_wc_rk
[result_TS_prob, alpha_prob] = KP_TSBGA(KP_uc_ak,1000,true,reps);
[result_TS_random_prob, ~] = KP_TSBGA(KP_uc_ak,1000,true,reps,true);
result_GA_prob = KP_TSBGA(KP_uc_ak,1000,false,reps);
result_CBR_prob = KP_CBRGA(KP_uc_ak,1000,reps);
result_TS{2} = result_TS_prob;
alpha{2} = alpha_prob;
result_TS_random{2} = result_TS_random_prob;
result_GA{2} = result_GA_prob;
result_CBR{2} = result_CBR_prob;

% KP_sc_rk
% [result_TS_prob, alpha_prob] = KP_TSBGA(KP_sc_ak,1000,true,reps);
% result_GA_prob = KP_TSBGA(KP_sc_ak,1000,false,reps);
% result_CBR_prob = KP_CBRGA(KP_sc_ak,1000,reps);
% result_TS{3} = result_TS_prob;
% alpha{3} = alpha_prob;
% result_GA{3} = result_GA_prob;
% result_CBR{3} = result_CBR_prob;

% result_DT = dependencyTree(KP_wc_ak,1000, reps);

% plot(mean(result_TS_prob))
% hold on 
% plot(mean(result_GA_prob))
% plot(mean(result_CBR_prob))
% % plot(mean(result_DT))
% hold off

% save('result_KP.mat', 'result_TS', 'alpha', 'result_GA', 'result_CBR');

%%
reps = 30;
for i = 1:length(result_TS)
  m_TS = mean(result_TS{i}); 
%   m_TS_random = mean(result_TS_random{i}); 
  m_GA = mean(result_GA{i}); 
  m_CBR = mean(result_CBR{i}); 

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

  box(axes1,'on');
  % Create legend
  legend1 = legend('AMTEA', 'CMA', 'TCIEA');
  set(legend1, 'FontSize', 15, 'Location', 'southeast', 'interpreter','latex');
  xlabel('Number of Function Evaluations', 'FontSize', 15, 'interpreter','latex')
  ylabel('Averaged Objective Value', 'FontSize', 15, 'interpreter','latex')
  
  if i == 1
    ylim([3000,4000])
  else
    ylim([2700,4200])
  end
  xlim([0, 5000])
  
  alpha_prob = alpha{i};
  alpha_Matrix = zeros(50,5);
  for rep = 1:reps
    alpha_Matrix = alpha_Matrix + alpha_prob{rep};
  end
  alpha_Matrix = alpha_Matrix/reps;
  
  error_Matrix = zeros(50, 5);
  for rep = 1:reps
      error_Matrix = error_Matrix + (alpha_prob{rep} - alpha_Matrix).^2;
  end
  error_Matrix = sqrt(error_Matrix/reps);

  figure2 = figure;

  % Create axes
  axes2 = axes('Parent',figure2);
  hold(axes2,'on');
  x = 2:2:100;
  plot2 = plot(x, alpha_Matrix(:,1),x, alpha_Matrix(:,2),...
      x, alpha_Matrix(:,3),x, alpha_Matrix(:,4),'Parent',axes2);
  set(plot2(1),'LineWidth',2,'Color',[0 0 1],'Marker','o', 'MarkerFaceColor',[0 0 1]);
  set(plot2(2),'LineWidth',2,'Color',[1 0 0],'Marker','v', 'MarkerFaceColor', [1 0 0]);
  set(plot2(3),'LineWidth',2,'Color',[0.5 0.5 0],'Marker','+', 'MarkerSize', 8);
  set(plot2(4),'LineWidth',2,'Color',[0 0 0],'Marker','square', 'MarkerFaceColor', [0 0 0]);
  shadedErrorBar(x, alpha_Matrix(:, 1)', error_Matrix(:, 1)','lineprops',{'Color',[0 0 1]},'transparent',1);
  shadedErrorBar(x, alpha_Matrix(:, 2)', error_Matrix(:, 2)','lineprops',{'Color',[1 0 0]},'transparent',1);
  shadedErrorBar(x, alpha_Matrix(:, 3)', error_Matrix(:, 3)','lineprops',{'Color',[0.5 0.5 0]},'transparent',1);
  shadedErrorBar(x, alpha_Matrix(:, 4)', error_Matrix(:, 4)','lineprops',{'Color',[0 0 0]},'transparent',1);

  ylim([0,0.8])
  legend2 = legend('KP\_uc\_rc (source 1)', 'KP\_wc\_rc (source 2)', 'KP\_sc\_rc (source 3)', 'KP\_sc\_ac (source 4)');
  set(legend2, 'FontSize', 15, 'interpreter','latex');
  xlabel('Number of Generations', 'FontSize', 15, 'interpreter','latex')
  ylabel('Transfer Coefficients ($\alpha_k$''s)', 'FontSize', 15, 'interpreter','latex')
end
