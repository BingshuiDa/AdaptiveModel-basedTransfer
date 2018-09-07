function [fitness_hist, alpha] = TSBGA(problem, dims, reps, trans, buildmodel, th_best) % Transfer stacking enhanced simple binary GA
% transfer = Boolean variable indicating if TRANSFER is to be used or not.
    load allmodels_TS
    if nargin < 5
        buildmodel = false;
        th_best = [];
    end
    if buildmodel == true
        load allmodels_CBR
    end
    
    pop = 200;
    gen = 1000;
    transfer = trans.transfer;
    if transfer
        TrInt = trans.TrInt; % Transfer interval ==> Generations after which transfer stacking is repeated
    end
    fitness_hist = zeros(reps, gen);
    bestsol = [];
    
    for rep = 1:reps
        alpha_rep = [];
        population = round(rand(pop,dims));
        fitness = funceval(population,problem,dims);
        [best_fit, ind] = max(fitness);
        disp(['Max fitness = ',num2str(best_fit)]); 
        fitness_hist(rep, 1) = max(fitness);
        bestsol = population(ind, :);
        
        for i = 2:gen
            if transfer && mod(i-1,TrInt) == 0
                mmodel = MixtureModel(allmodels_TS);
                mmodel = MixtureModel.createtable(mmodel,population,true,'umd');
                mmodel = MixtureModel.EMstacking(mmodel); %Recombination of probability models
                mmodel = MixtureModel.mutate(mmodel); % Mutation of stacked probability model
                offspring = MixtureModel.sample(mmodel,pop);
                alpha_rep = [alpha_rep; mmodel.alpha];
%             elseif mod(i-1,TrInt) == 0 && ~transfer
%                 allmodels_TS = {};
%                 mmodel = MixtureModel(allmodels_TS);
%                 mmodel = MixtureModel.createtable(mmodel,population,true,'umd');
%                 mmodel = MixtureModel.EMstacking(mmodel); %Recombination of probability models
%                 mmodel = MixtureModel.mutate(mmodel); % Mutation of stacked probability model
%                 offspring = MixtureModel.sample(mmodel,pop);
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
            cfitness = funceval(offspring,problem,dims);
            interpop = [population;offspring];
            interfitness = [fitness;cfitness];
            [interfitness,index] = sort(interfitness,'descend');
            fitness = interfitness(1:pop);
            interpop = interpop(index,:);        
            population = interpop(1:pop,:);
%             % RouletteWheel
%             for k = 1:pop
%                 index = RouletteWheelSelection(interfitness);
%                 population(k,:) = interpop(index,:);
%                 fitness(k) = interfitness(index);
%             end
%             % Binary tounament
%             ind = randperm(2*pop);
%             for k = 1:pop
%                 if interfitness(ind(2*k - 1)) > interfitness(ind(2*k))
%                     population(k,:) = interpop(ind(2*k - 1),:);
%                     fitness(k) = interfitness(ind(2*k - 1));
%                 else
%                     population(k,:) = interpop(ind(2*k),:);
%                     fitness(k) = interfitness(ind(2*k));
%                 end
%             end
            disp(['Max fitness = ',num2str(max(fitness))]); 
            fitness_hist(rep, i) = max(fitness);
            if buildmodel && sum(bestsol == population(1,:))~= dims
                allmodels_CBR = [allmodels_CBR; population(1,:)];
            end
            bestsol = population(1,:);
            
            if  buildmodel && fitness(1) == th_best
                model = ProbabilityModel('umd');
                noise = [];
                noise = round(rand(0.1*pop,dims));
                model = ProbabilityModel.buildmodel(model,[population;noise]);
                allmodels_TS{length(allmodels_TS)+1} = model;
                save('allmodels_TS.mat','allmodels_TS')
                buildmodel = false;

                allmodels_CBR = [allmodels_CBR; population(1,:)];
                break;
            end
        end 
        alpha{rep} = alpha_rep;
    end
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