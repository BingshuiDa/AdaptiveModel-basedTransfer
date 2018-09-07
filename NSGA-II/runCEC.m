%% Model preparation
clear 
clc

allmodels_TS = {}; % transfer stacking models
save('allmodels_TS','allmodels_TS')

allmodels_CBR = []; % Case-Based reasoning
save('allmodels_CBR','allmodels_CBR')

% UF1
name = 'UF1';
f = cec09(name);
dim = 30;
xrange = xboundary(name, dim);
L = xrange(:,1)';
U = xrange(:,2)';
pop = 50;
no_of_objs = 2;
NSGA_II(f, no_of_objs, L, U, pop);

% UF4
name = 'UF4';
f = cec09(name);
dim = 30;
xrange = xboundary(name, dim);
L = xrange(:,1)';
U = xrange(:,2)';
pop = 50;
no_of_objs = 2;
NSGA_II(f, no_of_objs, L, U, pop);


%% Transfer Stacking
close all
name = 'UF7';
f = cec09(name);
dim = 30;
xrange = xboundary(name, dim);
L = xrange(:,1)';
U = xrange(:,2)';
pop = 50;
no_of_objs = 2;

result_TS = [];
alpha = [];
result_NSGA = [];
result_CBR = [];
for rep = 1:30
    disp('Transfer stacking')
    [score_TS, alpha_TS] = TS_NSGA(f, no_of_objs, L, U, pop, true);
    result_TS = [result_TS; score_TS];
    alpha = [alpha; alpha_TS];

    disp('No transfer')
    score_NSGA = TS_NSGA(f, no_of_objs, L, U, pop, false);
    result_NSGA = [result_NSGA; score_NSGA];

    disp('Case-injected')
    score_CBR = CBR_NSGA(f, no_of_objs, L, U, pop);
    result_CBR = [result_CBR; score_CBR];
end