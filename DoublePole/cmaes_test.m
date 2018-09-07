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

% settings for cmaes
% opts.LBounds = -6;
% opts.UBounds = 6;
opts.StopFunEvals = 10000;
opts.StopFitness = -1999.99;
% opts.PopSize = 100;
% result_CMAES = [];
for rep = 1:5
    xstart = 12*rand(dim, 1) - 6;
    
    [~, ~, ~, ~, out, ~] = cmaes(f, xstart, 4, opts);
    disp([out.bestever, out.evals, out.stopflag])
end