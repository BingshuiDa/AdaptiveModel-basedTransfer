function model = BGA(problem,dims,th_best) % Simple Binary GA, with uniform crossover and bit-flip mutation
% Function returns a probability distribution model encapsulating the structure of the
% favorable regions of the search space.
    load allmodels_TS % allmodels should be initialized as an empty cell array where source probability distributions are gradually stored
%     load allmodels_CBR

    pop = 200;
    gen = 1000;
    
    population = round(rand(pop,dims));    
    fitness = funceval(population,problem,dims);
    buildmodel = true;
    bestfitness = max(fitness);
    
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
        cfitness = funceval(offspring,problem,dims);
        interpop = [population;offspring];
        interfitness = [fitness;cfitness];
        [interfitness,index] = sort(interfitness,'descend');
        fitness = interfitness(1:pop);
        interpop = interpop(index,:);        
        population = interpop(1:pop,:);
        disp(['Best fitness = ',num2str(fitness(1))]);
%         if bestfitness ~= fitness(1)
%             allmodels_CBR = [allmodels_CBR; population(1,:)];
%         end
        bestfitness = fitness(1);
        if (fitness(1) >= th_best || i == gen) && buildmodel
            model = ProbabilityModel('umd');
%             noise = round(rand(0.1*pop,dims));
            model = ProbabilityModel.buildmodel(model,population);
            allmodels_TS{length(allmodels_TS)+1} = model;
            save('allmodels_TS.mat','allmodels_TS')
            buildmodel = false;
            
%             allmodels_CBR = [allmodels_CBR; population(1,:)];
            break;
        end
    end
%     save('allmodels_CBR.mat', 'allmodels_CBR');
end

function fitness = funceval(population,problem,dims)
    if strcmp(problem,'onemax')
        fitness = sum(population,2);
    elseif strcmp(problem,'onemin')
        fitness = dims - sum(population,2);
    elseif strcmp(problem,'trap5')
        pop = size(population, 1);
        fitness = zeros(pop,1);
        index = 1:dims;
        index = vec2mat(index,5);
        rows = size(index,1);
        for i = 1:pop
            fitsum = 0;
            for j = 1:rows
                contri = sum(population(i,index(j,:)));
                if contri == 5
                    fitsum = fitsum+5;
                else
                    fitsum = fitsum+(4-contri);
                end
            end
            fitness(i) = fitsum;
        end
    elseif strcmp(problem,'trap2')
        pop = size(population, 1);
        fitness = zeros(pop,1);
        index = 1:dims;
        index = vec2mat(index,2);
        rows = size(index,1);
        for i = 1:pop
            fitsum = 0;
            for j = 1:rows
                contri = sum(population(i,index(j,:)));
                if contri == 2
                    fitsum = fitsum+2;
                else
                    fitsum = fitsum+(1-contri);
                end
            end
            fitness(i) = fitsum;
        end
    elseif strcmp(problem,'hamming')
        load 'hamming_optimal.mat'
        optimal = optimal(1:dims);
        pop = size(population, 1);
        fitness = dims - sum(abs(population - repmat(optimal, pop, 1)), 2);
    end
end