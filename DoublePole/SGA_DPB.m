function [bestSol, bestfit_hist] = SGA_DPB(h_len)
load allmodels_TS
load allmodels_CBR

pop = 100;
generation = 100;
bestSol = [];

muc = 20;
mum = 20;

bestfit_hist = zeros(1,generation);

% initiailize cart 
cart = Cart();

% initilize network
nInput = 6; nHidden = 10; nOutput = 1;
net.size(1) = nInput;
net.size(2) = nHidden;
net.size(3) = nOutput;
net.transferFcn{1} = 'tansig';
net.transferFcn{2} = 'tansig';

nHidVariables = nInput*nHidden;
nOutVariables = nHidden*nOutput;
nVariable = nHidVariables + nOutVariables;

population = 12*rand(pop, nVariable) - 6;
fitness = zeros(1,pop);
for i = 1:pop
  net.IW = reshape(population(i,1:nHidVariables), nHidden, nInput);
  net.LW = population(i,(nHidVariables+1):end);
  fitness(i) = evaluate(net, cart, h_len);
end

bestfit_hist(1) = max(fitness);
disp(['Generation 1: ', num2str(bestfit_hist(1))]);

for gen = 2:generation
  randlist=randperm(pop);
  child = zeros(size(population));
  for i = 1:2:pop-1
    p1=randlist(i);
    p2=randlist(i+1);
    [child(i,:), child(i+1,:)]=crossover(population(p1,:), population(p2,:),muc,nVariable);
    child(i,:) = mutate(child(i,:),mum,nVariable);
    child(i+1,:)=mutate(child(i+1,:),mum,nVariable);
  end
  
  cfitness = zeros(1,pop);
  for i = 1:pop
    net.IW = reshape(child(i,1:nHidVariables), nHidden, nInput);
    net.LW = child(i,(nHidVariables+1):end);
    cfitness(i) = evaluate(net, cart, h_len);
  end  
  
  intpopulation = [population; child];
  intfitness = [fitness, cfitness];
  [bestfit_hist(gen),ind] = max(intfitness);
  bestSol = [bestSol; intpopulation(ind,:)];
  
  % elitest
  [~, ind] = sort(intfitness, 2, 'descend');
  fitness = intfitness(ind(1:pop));
  population = intpopulation(ind(1:pop), :);
  
  
  % binary tournament
%   randlist = randperm(2*pop);
%   count = 1;
%   for i = 1:2:2*pop-1
%     if intfitness(randlist(i)) > intfitness(randlist(i+1))
%       population(count,:) = intpopulation(randlist(i),:);
%       fitness(count) = intfitness(randlist(i));
%     else
%       population(count,:) = intpopulation(randlist(i+1),:);
%       fitness(count) = intfitness(randlist(i+1));
%     end
%     count = count + 1;
%   end
  disp(['Generation ', num2str(gen), ': ',num2str(bestfit_hist(gen))]);
  
  if (bestfit_hist(gen) - 2000) > -0.0001
    disp('Found solution!')
    bestfit_hist((gen+1):end) = bestfit_hist(gen);
    break;
  end
end
if bestfit_hist(end) - 2000 > -0.0001
  model  = ProbabilityModel('mvarnorm');
  noisypop = [];
  for i = 1:1.0*pop 
    if i<= pop
      noisypop = [noisypop; population(i,:)]; 
    else
      noisypop = [noisypop; rand(1, nVariable)*12-6];
    end
  end
  model = ProbabilityModel.buildmodel(model,noisypop);
  allmodels_TS{length(allmodels_TS)+1} = model;
  save('allmodels_TS.mat', 'allmodels_TS')
  allmodels_CBR = [allmodels_CBR; bestSol];
  save('allmodels_CBR.mat', 'allmodels_CBR')
end
end

function time = evaluate(net, cart, h_len)
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
end

function val = eval_net(net, input)
val_Hid = feval(net.transferFcn{1}, net.IW*input);
val = feval(net.transferFcn{2}, net.LW*val_Hid);
end

function [rnvec1,rnvec2]=crossover(p1,p2,muc,dim)
rnvec1=zeros(1,dim);
rnvec2=zeros(1,dim);
randlist=rand(1,dim);
for i=1:dim
    if randlist(i)<=0.5
        k=(2*randlist(i))^(1/(muc+1));
    else
        k=(1/(2*(1-randlist(i))))^(1/(muc+1));
    end
    rnvec1(i)=0.5*(((1 + k)*p1(i)) + (1 - k)*p2(i));
    rnvec2(i)=0.5*(((1 - k)*p1(i)) + (1 + k)*p2(i));
    if rand(1) < 0.5
      tmp = rnvec1(i);
      rnvec1(i) = rnvec2(i);
      rnvec2(i) = tmp;
    end
end
end

function rnvec=mutate(p,mum,dim)
rnvec=p;
for i=1:dim
    if rand(1)<1/dim
        u=rand(1);
        if u <= 0.5
            del=(2*u)^(1/(1+mum)) - 1;
            rnvec(i)=p(i) + del*(p(i));
        else
            del= 1 - (2*(1-u))^(1/(1+mum));
            rnvec(i)=p(i) + del*(1-p(i));
        end
    end
end

% for i=1:dim
%     if rand(1)<0.2
%         rnvec(i) = rand(1)*12-6;
%     end
% end
end