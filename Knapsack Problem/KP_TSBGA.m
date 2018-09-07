function [fitness_hist, alpha] = KP_TSBGA(problem,dims,transfer,reps,random_coef) % Transfer stacking enhanced simple binary GA
% transfer = Boolean variable indicating if TRANSFER is to be used or not.
    if nargin < 5
      random_coef = false;
    end
    load allmodels_TS
%     allmodels{length(allmodels)+1} = ProbabilityModel('umd');

    pop = 50;
    gen = 100;
    TrInt = 2; % Transfer interval ==> Generations after which transfer stacking is repeated
    alpha = {};
    fitness_hist = zeros(reps, gen);
    
    for rep = 1:reps
        population = round(rand(pop,dims));
        fitness = funceval(population,problem,dims,pop);
        alpha_rep = [];
        
        for i = 1:gen
            if mod(i,TrInt) == 0 && transfer
    % %             model = ProbabilityModel('umd');
    % %             noise = round(rand(round(0.1*pop),dims)); % Random noise term to prevent overfitting to target optimiztion task
    % %             model = ProbabilityModel.buildmodel(model,[population;noise]);
    % %             allmodels{length(allmodels)} = model;
                mmodel = MixtureModel(allmodels_TS);
                mmodel = MixtureModel.createtable(mmodel,population,true,'umd');
                mmodel = MixtureModel.EMstacking(mmodel, random_coef); %Recombination of probability models
                mmodel = MixtureModel.mutate(mmodel); % Mutation of stacked probability model
                offspring = MixtureModel.sample(mmodel,pop);
                alpha_rep = [alpha_rep; mmodel.alpha];
            else
                parent1 = population(randperm(pop),:);
                parent2 = population(randperm(pop),:);
                tmp = rand(pop,dims);        
                offspring = zeros(pop,dims);
                index = tmp>=0.5;
                offspring(index) = parent1(index);
                index = tmp<0.5;
                offspring(index) = parent2(index);
                tmp = rand(pop,dims);
                index = tmp<(1/dims);
                offspring(index) = abs(1-offspring(index));
            end        
            cfitness = funceval(offspring,problem, dims, pop);
            interpop = [population;offspring];
            interfitness = [fitness;cfitness];
            [interfitness,index] = sort(interfitness,'descend');
            fitness = interfitness(1:pop);
            interpop = interpop(index,:);        
            population = interpop(1:pop,:);
            disp(['Generation = ' num2str(i) ' ' 'Mean fitness = ',num2str(mean(fitness))]); 
            fitness_hist(rep, i) = mean(fitness);
%             if fitness(1) == problem.opt && premSTOP
%                 break;
%             end
        end  
        alpha{rep} = alpha_rep;
    end
%     population(1,:)
end

function fitness = funceval(population,problem,dims,pop)
    fitness = zeros(pop,1);
    Weights = problem.w';
    Profits = problem.p';
    Ratios =    Profits./Weights;
    for i = 1:pop
        BV = population(i,:);
        TotalWeight = sum(BV*Weights);
        TotalProfit = sum(BV*Profits);
        
        if TotalWeight > problem.cap % Repair solution
            selections = sum(BV);
            List = zeros(selections,2);
            counter = 1;
            for j = 1 : dims
                if BV(j) == 1
                    List(counter,1) = Ratios(j);
                    List(counter,2) = j;
                    counter = counter + 1;
                end
                if counter > selections
                    break;
                end
            end
            List = sortrows(List,-1);
            counter = selections;
            while TotalWeight > problem.cap
                BV(List(counter,2)) = 0;
                TotalWeight = TotalWeight - Weights(List(counter,2));
                TotalProfit = TotalProfit - Profits(List(counter,2));
                counter = counter - 1;
            end
        end
        fitness(i) = TotalProfit;
    end
end