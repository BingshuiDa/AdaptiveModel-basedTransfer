function fitness_hist = CBRBGA(problem,dims,reps,TrInt) % Transfer stacking enhanced simple binary GA
% transfer = Boolean variable indicating if TRANSFER is to be used or not.
    load allmodels_CBR
%     allmodels = [];

    pop = 50;
    injection = 0.1*pop;
    gen = 100;
    TrInt = 2; % Transfer interval ==> Generations after which transfer stacking is repeated
    fitness_hist = zeros(reps, gen);
    
    
    for rep = 1:reps
        population = round(rand(pop,dims));
%         randlist = randperm(size(allmodels_CBR,1));
%         population(pop,:) = allmodels_CBR(randlist(1),:);
        fitness = funceval(population,problem,dims);
        disp(['Max fitness = ',num2str(max(fitness))]); 
        fitness_hist(rep, 1) = max(fitness);

        for i = 2:gen
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
            
            if mod(i-1,TrInt) == 0
                dist = sq_dist(population(1,:)', allmodels_CBR');
                [~, ind] = sort(dist);
                offspring((pop-injection+1):pop,:) = allmodels_CBR(ind(1:injection),:);
            end
            
            cfitness = funceval(offspring,problem,dims);
            interpop = [population;offspring];
            interfitness = [fitness;cfitness];
            [interfitness,index] = sort(interfitness,'descend');
            fitness = interfitness(1:pop);
            interpop = interpop(index,:);        
            population = interpop(1:pop,:);
            disp(['Mean fitness = ',num2str(max(fitness))]); 
            fitness_hist(rep, i) = max(fitness);
        end 
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
    end        
end