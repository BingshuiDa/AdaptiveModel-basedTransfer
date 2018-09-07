allmodels_TS = {};
save('allmodels_TS.mat', 'allmodels_TS')
allmodels_CBR = [];
save('allmodels_CBR.mat', 'allmodels_CBR')

len = [0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.775];

for i = 1:length(len)
    cart = Cart();
    % half length
    h_len = len(i)/2;
    % define network
    nInput = 6; nHidden = 10; nOutput = 1;
    net.size(1) = nInput;
    net.size(2) = nHidden;
    net.size(3) = nOutput;
    net.transferFcn{1} = 'tansig';
    net.transferFcn{2} = 'tansig';
    dim = nInput*nHidden + nHidden*nOutput;

    f = @(x)eval_cart(x, net, cart, h_len);

    while 1
        x = 12*rand(1, dim) - 6;
        if i > 11
            
            load allmodels_CBR.mat
            x = allmodels_CBR(end, :);
        end
        [xopt, fopt, neval] = xnes_source(f, dim, x', 10000);
        disp([fopt, neval])
        if fopt < (-2000 + 0.0001)
            disp(['Source task for length ', num2str(len(i)), ' successfully built!'])
            break;
        end
    end
end