function val = eval_net(net, input)
val_Hid = feval(net.transferFcn{1}, net.IW*input);
val = feval(net.transferFcn{2}, net.LW*val_Hid);
end