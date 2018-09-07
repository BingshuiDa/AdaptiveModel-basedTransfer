function model = KP_BGA(problem,dims,premSTOP) % Simple Binary GA, with uniform crossover and bit-flip mutation for the 0/1 Knapsack problem
% Function returns a probability distribution model encapsulating the structure of the
% favorable regions of the search space.
    load allmodels_TS; % allmodels should be initialized as an empty cell array where source probability distributions are gradually stored
    load allmodels_CBR
    
    pop = 100;
    gen = 100;
    
    population = round(rand(pop,dims));    
    fitness = funceval(population,problem,dims,pop);
    bestfitness = fitness(1);
    counter = 0;
    for i = 1:gen
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
        cfitness = funceval(offspring,problem,dims,pop);
        interpop = [population;offspring];
        interfitness = [fitness;cfitness];
        [interfitness,index] = sort(interfitness,'descend');
        fitness = interfitness(1:pop);
        interpop = interpop(index,:);        
        population = interpop(1:pop,:);
        disp(['Best fitness = ',num2str(fitness(1))]);  
        allmodels_CBR = [allmodels_CBR; population(1,:)]; % save in case base
        if fitness(1) > bestfitness
            bestfitness = fitness(1);
            counter = 0;
        else
            counter = counter+1;
        end
        if counter == 20 && premSTOP
            break;
        end
    end
    
    allmodels_CBR = [allmodels_CBR; population(1,:)]; % save in case base
    save('allmodels_CBR.mat', 'allmodels_CBR');
    
    problem.opt
    model = ProbabilityModel('umd');
    noise = [];
    noise = round(rand(0.1*pop,dims));
    model = ProbabilityModel.buildmodel(model,[population;noise]);
    allmodels_TS{length(allmodels_TS)+1} = model;
    save('allmodels_TS','allmodels_TS')
end

function fitness = funceval(population,problem,dims,pop)
    fitness = zeros(pop,1);
    Weights = problem.w';
    Profits = problem.p';
    Ratios = Profits./Weights;
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