result_TS = [];
result_GA = [];
result_CBR = [];
result_NES = [];
neval_NES = [];
reps = 50;
h_len = 0.4;

cart = Cart();
% define network
nInput = 6; nHidden = 10; nOutput = 1;
net.size(1) = nInput;
net.size(2) = nHidden;
net.size(3) = nOutput;
net.transferFcn{1} = 'tansig';
net.transferFcn{2} = 'tansig';
dim = nInput*nHidden + nHidden*nOutput;

f = @(x)eval_cart(x, net, cart, h_len);

for rep = 1:reps
  disp(rep);
  
  disp('Transfer Stacking')
  [~, bestfit_hist, alpha_rep] = TS_DPB(true, h_len);
  result_TS = [result_TS; bestfit_hist]; % transfer stacking
  alpha{rep} = alpha_rep;
  
  disp('GA')
  [~, bestfit_hist] = TS_DPB(false, h_len);
  result_GA = [result_GA; bestfit_hist]; % naive GA
  
  disp('CBR')
  [~, bestfit_hist] = CBR_DPB(h_len); % case-injected GA
  result_CBR = [result_CBR; bestfit_hist];

  disp('NES')
  x = 12*rand(1, dim) - 6;
  [xopt, fopt, neval] = xnes(f, dim, x', 10000);
  disp([fopt, neval])
  result_NES = [result_NES, fopt];
  neval_NES = [neval_NES, neval];
  
%   disp('CMAES')
%   opts.StopFunEvals = 1000;
%   opts.PopSize = 20;
%   opts.StopFitness = -1999.9;
%   nInput = 6; nHidden = 10; nOutput = 1;
%   nHidVariables = nInput*nHidden;
%   nOutVariables = nHidden*nOutput;
%   nVariable = nHidVariables + nOutVariables;
%   xstart = rand(nVariable,opts.PopSize).*12 + 6;
%   [~,~,~,~,out,~] = cmaes('fnceval', xstart, [], opts);
%   result_CMAES = [result_CMAES; out.meanfitness];
end

% save('result_DP5.mat', 'result_TS', 'alpha', 'result_GA', 'result_CBR', 'result_NES');

%%
YMatrix = [mean(result_TS); mean(result_GA); mean(result_CBR)]';
figure1 = figure1;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot1 = plot(YMatrix,'LineWidth',1,'Parent',axes1);
set(plot1(1),'DisplayName','GA+TO','Color',[0 0 1]);
set(plot1(2),'DisplayName','GA','Color',[1 0 0]);
set(plot1(3),'DisplayName','CBR','Color',[0 0 0]);
box(axes1,'on');
% Create legend
legend(axes1,'show');

%%
success_TS = find(result_TS(:,end) >= (2000 - 0.001));
success_GA = find(result_GA(:,end) >= (2000 - 0.001));
success_CBR = find(result_CBR(:,end) >= (2000 - 0.001));

nfe_TS = [];
for i = 1:length(success_TS)
  nfe_TS = [nfe_TS, min(find(result_TS(success_TS(i),:) >= (2000 - 0.001))) * 100];
end

nfe_GA = [];
for i = 1:length(success_GA)
  nfe_GA = [nfe_GA, min(find(result_GA(success_GA(i),:) >= (2000 - 0.001))) * 100];
end

nfe_CBR = [];
for i = 1:length(success_CBR)
  nfe_CBR = [nfe_CBR, min(find(result_CBR(success_CBR(i),:) >= (2000 - 0.001))) * 100];
end


%%
success_TS = find(result_TS(:,end) >= (2000 - 0.001));
mean_alpha = zeros(5,10);
std_alpha = zeros(5,10);

alpha_sample = {};
for i = 1:10
    alpha_sample{i} = [];
    for j = 1:length(success_TS)
        ind = success_TS(j);
        if size(alpha{ind}, 1) < i
            continue;
        end
        alpha_sample{i} = [alpha_sample{i}, alpha{ind}(i,9)];
    end
end
for i = 1:10
    mean_alpha(1,i) = mean(alpha_sample{i});
    std_alpha(1,i) = std(alpha_sample{i});
end

alpha_sample = {};
for i = 1:10
    alpha_sample{i} = [];
    for j = 1:length(success_TS)
        ind = success_TS(j);
        if size(alpha{ind}, 1) < i
            continue;
        end
        alpha_sample{i} = [alpha_sample{i}, alpha{ind}(i,10)];
    end
end
for i = 1:10
    mean_alpha(2,i) = mean(alpha_sample{i});
    std_alpha(2,i) = std(alpha_sample{i});
end

alpha_sample = {};
for i = 1:10
    alpha_sample{i} = [];
    for j = 1:length(success_TS)
        ind = success_TS(j);
        if size(alpha{ind}, 1) < i
            continue;
        end
        alpha_sample{i} = [alpha_sample{i}, alpha{ind}(i,11)];
    end
end
for i = 1:10
    mean_alpha(3,i) = mean(alpha_sample{i});
    std_alpha(3,i) = std(alpha_sample{i});
end

alpha_sample = {};
for i = 1:10
    alpha_sample{i} = [];
    for j = 1:length(success_TS)
        ind = success_TS(j);
        if size(alpha{ind}, 1) < i
            continue;
        end
        alpha_sample{i} = [alpha_sample{i}, alpha{ind}(i,12)];
    end
end
for i = 1:10
    mean_alpha(4,i) = mean(alpha_sample{i});
    std_alpha(4,i) = std(alpha_sample{i});
end

alpha_sample = {};
for i = 1:10
    alpha_sample{i} = [];
    for j = 1:length(success_TS)
        ind = success_TS(j);
        if size(alpha{ind}, 1) < i
            continue;
        end
        alpha_sample{i} = [alpha_sample{i}, alpha{ind}(i,13)];
    end
end
for i = 1:10
    mean_alpha(5,i) = mean(alpha_sample{i});
    std_alpha(5,i) = std(alpha_sample{i});
end

x = 2:10:100;

figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot1 = plot(x, mean_alpha(1,:),x, mean_alpha(2,:),...
    x, mean_alpha(3,:),x, mean_alpha(4,:),...
    x, mean_alpha(5,:),'Parent',axes1);
set(plot1(1),'LineWidth',2,'Color',[0 0 1],'Marker','o', 'MarkerFaceColor',[0 0 1]);
set(plot1(2),'LineWidth',2,'Color',[1 0 0],'Marker','v', 'MarkerFaceColor', [1 0 0]);
set(plot1(3),'LineWidth',2,'Color',[0.5 0.5 0],'Marker','+', 'MarkerSize', 8);
set(plot1(4),'LineWidth',2,'Color',[0 0 0],'Marker','square', 'MarkerFaceColor', [0 0 0]);
set(plot1(5),'LineWidth',2,'Color',[0.5 0 0.5],'Marker','x', 'MarkerSize', 8);
shadedErrorBar(x, mean_alpha(1, :), 0.2*std_alpha(1, :),'lineprops',{'Color',[0 0 1]},'transparent',1);
shadedErrorBar(x, mean_alpha(2, :), 0.2*std_alpha(2, :),'lineprops',{'Color',[1 0 0]},'transparent',1);
shadedErrorBar(x, mean_alpha(3, :), 0.2*std_alpha(3, :),'lineprops',{'Color',[0.5 0.5 0]},'transparent',1);
shadedErrorBar(x, mean_alpha(4, :), 0.2*std_alpha(4, :),'lineprops',{'Color',[0 0 0]},'transparent',1);
shadedErrorBar(x, mean_alpha(5, :), 0.2*std_alpha(5, :),'lineprops',{'Color',[0.5 0 0.5]},'transparent',1);

ylim([0,0.35])
box(axes1,'on');
legend1 = legend( '$e_{0.6}$ (source 1)', '$e_{0.65}$ (source 2)', '$e_{0.7}$ (source 3)', ...
    '$e_{0.75}$ (source 4)', '$e_{0.775}$ (source 5)');
set(legend1, 'FontSize', 15, 'interpreter','latex');
xlabel('Number of Generations', 'FontSize', 15, 'interpreter','latex')
ylabel('Transfer Coefficient ($\alpha_k$''s)', 'FontSize', 15, 'interpreter','latex')


