function val = eval_cart(input, net, cart, h_len)
input = reshape(input, 1, numel(input));
nInput = net.size(1);
nHidden = net.size(2);
nOutput = net.size(3);

nHidVar = nInput*nHidden;

net.IW = reshape(input(1:nHidVar), nHidden, nInput);
net.LW = input((nHidVar+1):end);

cart = initialize(cart, h_len);

while 1
  state = get_state(cart);
  cart.applied_force = 10*eval_net(net, state);
  cart = update_state(cart);
  cart = update_state(cart);
  if cart.failed 
    time = cart.time;
    break;
  elseif (cart.time - 2000) > -0.00001
      time = cart.time;
      break;
  end
end
val = -time;
end

function val = eval_net(net, input)
val_Hid = feval(net.transferFcn{1}, net.IW*input);
val = feval(net.transferFcn{2}, net.LW*val_Hid);
end