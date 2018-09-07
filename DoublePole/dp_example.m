neval_NES = [];
reps = 1;
h_len = 0.1;

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
	x = 12*rand(1, dim) - 6;
	[xopt, fopt, neval] = xnes(f, dim, x', 10000);
	disp([fopt, neval])
	result_NES = [result_NES, fopt];
	neval_NES = [neval_NES, neval];
end
