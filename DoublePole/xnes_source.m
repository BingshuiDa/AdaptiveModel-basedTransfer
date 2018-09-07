function [xopt,fopt,neval] = xnes_source(f,d,x,num_eval)
% Written by Sun Yi (yi@idsia.ch).
% Modified by Bingshui Da

% parameters
% L = 4+3*floor(log(d));
load allmodels_TS
load allmodels_CBR

bestSol = [];
L = 100;
etax = 1; etaA = 0.5*min(1.0/d,0.25);
shape = max(0.0, log(L/2+1.0)-log(1:L)); shape = shape / sum(shape);

% initialize
xopt = x; 
fopt = f(x);
neval = 1;
A = zeros(d);
weights = zeros(1,L);
fit = zeros(1,L);
% tm = cputime;
gen = 1;

while neval < num_eval
    expA = expm(A);
    
    % step 1: sampling & importance mixing
    Z = randn(d,L); X = repmat(x,1,L)+expA*Z;
    for i = 1 : L, fit(i) = f(X(:,i)); end
    neval = neval + L;
    
    % step 2: fitness reshaping
    [~, idx] = sort(fit); weights(idx) = shape;
    if fit(idx(1)) < fopt
        xopt = X(:,idx(1)); fopt = fit(idx(1));
        bestSol = [bestSol; xopt'];
    end
    disp(['Generation: ', num2str(gen), ', current best: ', num2str(fopt)]);

    % step 3: compute the gradient for C and x
    G = (repmat(weights,d,1).*Z)*Z' - sum(weights)*eye(d);
    dx = etax * expA * (Z*weights');
    dA = etaA * G;
  
    % step 4: compute the update  
    x = x + dx; A = A + dA;
    gen = gen + 1;
    
    if trace(A)/d < -10*log(10), break; end
    if fopt < (-2000 + 0.0001)
        disp(['Solution found!   ', num2str(neval)]);
        % save CBR
        allmodels_CBR = [allmodels_CBR; bestSol];
        save('allmodels_CBR.mat', 'allmodels_CBR')
        % build probability source model
        Z = randn(d,L); X = repmat(x,1,L)+expA*Z;
        model = ProbabilityModel('mvarnorm');
        model = ProbabilityModel.buildmodel(model, X');
        allmodels_TS{length(allmodels_TS)+1} = model;
        save('allmodels_TS.mat', 'allmodels_TS')
        break;
    end
end
