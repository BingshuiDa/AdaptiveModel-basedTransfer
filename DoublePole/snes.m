function[mean_vec, var_vec, mean_vec_fit] = snes(fitness, num_iter, dim, pop_size, learn_rates)
%[mean_vec, var_vec, mean_vec_fit] = 
%           snes(dim, mean_vec, var_vec, pop_size, learn_rates, @fitness)
%INPUTS:
% @fitness: fitness function (e.g., @rosenfit for rosenfit.m)
%           which takes an individual as input and outputs the fitness
% num_iter: number of iterations
% dim: dimension
% **pop_size: number of units to evaluate 
% **learn_rates: vector of two, for mean and variance
% defaults if pop_size and learn_rates not supplied
%OUTPUTS
% mean_vec: final mean
% var_vec: final variance
% mean_vec_fit: the fitness of the mean, measured at each iteration
%
%code by Matt Luciw (matt.luciw at gmail)
%contact Tom Schaul with all your SNES questions

%initial mean and variance
mean_vec = rand(dim,1);
var_vec = ones(dim,1);  %ones!  <-- this is important

if (nargin < 4)
    %abra cadabra
    pop_size = 4 + floor(3 * log(dim));
    
    learn_rates = [1 (3 + log(dim))/(5 * sqrt(dim))];
end

mean_vec_fit = zeros(1,num_iter);

%outer loop: number of population evaluations
for i = 1:num_iter
     
    if (mod(i,500)==0) 
        fprintf(1, '\nGeneration %d...', i)
    end
    
    %draw from standard normal distribution
    curr_samples = randn(pop_size,dim);
    
    %add the input mean and variance
    curr_members = (curr_samples .* repmat(var_vec',pop_size,1)) + ...
        repmat(mean_vec',pop_size,1);
    
    %store samples
    S = curr_samples';
        
    %inner loop: number of population members
    for j = 1 : pop_size
        
        %fitness evaluated here for this sample (and stored)
        fit(j) = fitness(curr_members(j,:));
        
    end
    
    %sort by fitness so most fit guys are last
    [dummy order] = sort(fit);
    
    %ordered set of samples
    S = S(:,order);
    
    %utilities which must sum to one
    %first half of the population has zero utility
    threshold = floor(pop_size / 2);
    step_size = 1 / threshold; 
    U = zeros(1,pop_size);
    U(end-threshold+1:end) = step_size:step_size:1;
    U = U ./ sum(U);
    
    %compute gradients
    %one for mean
    mean_grad = U*S';

    %variance gradient
    S_sq_minus = S.^2 - 1;
    var_grad = U * S_sq_minus';

    %update parameters
    mean_vec = mean_vec + learn_rates(1) * var_vec .* mean_grad';
 
    var_vec = var_vec .* exp(learn_rates(2) / 2 * var_grad)';
    
    %evaluate fitness of mean (for plotting)
    mean_vec_fit(i) = fitness(mean_vec);
    
    %uncomment for spiffy updating plot
    %plot(mean_vec_fit)
    %drawnow
end


