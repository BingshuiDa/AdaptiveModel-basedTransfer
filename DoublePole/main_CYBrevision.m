%This is additional experiment for the CYB reviews. This is to compare with
%Bayesian optimization algorithm on double pole balancing problem uses the
%bayesopt in matlab Statistics and Machine learning Toolbox.

%% define objective function
addpath(genpath('bayesopt'))
cart = Cart();
h_len = 0.05;
% define network
nInput = 6; nHidden = 10; nOutput = 1;
net.size(1) = nInput;
net.size(2) = nHidden;
net.size(3) = nOutput;
net.transferFcn{1} = 'tansig';
net.transferFcn{2} = 'tansig';
dim = nInput*nHidden + nHidden*nOutput;

f = @(x)eval_cart(x, net, cart, h_len);

opt = defaultopt();
opt.dims = dim;
opt.mins = -6*ones(1,dim);
opt.maxes = 6*ones(1,dim);
opt.max_iters = 10000;
opt.optimize_ei = 1;
opt.grid_size = 300; % If you use the optimize_ei option
opt.do_cbo = 0;
opt.initial_points = 12*rand(1, dim) - 6;

fprintf('Optimizing hyperparamters of function "samplef.m" ...\n');
[ms,mv,T] = bayesopt(f,opt);   % ms - Best parameter setting found
                               % mv - best function value for that setting L(ms)
                               % T  - Trace of all settings tried, their function values, and constraint values.
                              
%% Print results
fprintf('******************************************************\n');
fprintf('Best hyperparameters:      P1=%2.4f, P2=%2.4f\n',ms(1),ms(2));
fprintf('Associated function value: F([P1,P2])=%2.4f\n',mv);
fprintf('******************************************************\n');

