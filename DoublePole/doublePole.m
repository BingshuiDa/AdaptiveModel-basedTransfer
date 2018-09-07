function time = doublePole(x, len)
%time=DOUBLEPOLE(x) : simulation of double pole balancing optimization
%problem as a black box. The structure of the NN is fixed, with no bias
%term in the network. The input length must be the same as the number of
%variables in the network. Maximum balancing time is set as 2000 seconds.
%
%INPUT:
%   x: weight for the feedforward neural network
%   len: the length of shorter pole ranging from 0.1 to 0.8
%       (the length of longer pole is fixed as 1. As len increases, the 
%   difficulty of the task increases. len=0.8 is extremely challenging.)
%
%OUTPUT:
%   time: the duration time the simulation can keep the double pole
%           balanced

h_len = len/2; 

cart = Cart();  % define the cart object
% define network, the network structure is fixed throughout
nInput = 6; nHidden = 10; nOutput = 1;
net.size(1) = nInput;
net.size(2) = nHidden;
net.size(3) = nOutput;
net.transferFcn{1} = 'tansig';
net.transferFcn{2} = 'tansig';
dim = nInput*nHidden + nHidden*nOutput; % the number of variables in the NN
if dim ~= length(x)
    error('Dimensionality mismatches!')
end
net.IW = reshape(x(1:nInput*nHidden), nHidden, nInput);
net.LW = reshape(x(nInput*nHidden + 1:end), 1, nHidden);
time = evaluate(net, cart, h_len); 
end

function time = evaluate(net, cart, h_len)
% evaluate the network to output the duration to balance the double pole
    cart = initialize(cart, h_len);
    while 1
      state = get_state(cart);
      cart.applied_force = 10*eval_net(net, state);
      cart = update_state(cart);
      cart = update_state(cart);
      if cart.failed || (cart.time - 2000) > -0.0001
        time = cart.time;
        break;
      end
    end
end

function val = eval_net(net, input)
% evaluate the network to output the force
    val_Hid = feval(net.transferFcn{1}, net.IW*input);
    val = feval(net.transferFcn{2}, net.LW*val_Hid);
end