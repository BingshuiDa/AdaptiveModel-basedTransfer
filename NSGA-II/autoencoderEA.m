function store = autoencoderEA(f, no_of_objs, L, U, pop, TrInt)
%store=AUTOENCODEREA(f, no_of_objs, L, U, pop)
% This is to implement autoencoder EA by Feng et al.
%
% Author: Bingshui Da <DA0002UI@e.ntu.edu.sg>

gen=100;
muc=20;% crossover index SBX
mum=2;% mutation index (polynomial mutation)

if mod(pop, 2) ~= 0
    pop = pop + 1;
end
dim = length(L);

store = [];
if ischar(f)
    if no_of_objs == 2
        data = load([f, '.pf']);
    elseif no_of_objs == 3
        data = load([f, '.3D.pf']);
    end
else
    data = load([f.name, '.pf']);
end

for i = 1:pop
    population(i) = Chromosome;
    population(i) = initialize(population(i), dim);
end

for i = 1:pop
    population(i) = evaluate(population(i), f, L, U);
end
[population, frontnumbers] = SolutionComparison.nondominatedsort(population, pop, no_of_objs);
[population, minimums] = SolutionComparison.diversity(population, frontnumbers, pop, no_of_objs);

for generation = 1:gen
    if  mod(generation, TrInt) == 0
        solutions = zeros(pop, dim);
        for i = 1:pop
            solutions(i,:) = population(i).rnvec;
        end
        [offspring, count] = denoiseAutoencode(solutions, generation);
        for i = 1:count
            child(i) = Chromosome;
            child(i).rnvec = offspring(i,:);
        end
        if count < pop
            % add child individuals
            if mod(pop - count, 2) ~= 0
                child_pop = pop - count + 1;
            else
                child_pop = pop - count;
            end
            for i = 1:child_pop
                parent(i) = Chromosome;
                p1 = 1+round(rand*(pop-1));
                p2 = 1+round(rand*(pop-1));
                if population(p1).front < population(p2).front
                    parent(i).rnvec = population(p1).rnvec;
                elseif population(p1).front == population(p2).front
                    if population(p1).CD > population(p2).CD
                        parent(i).rnvec = population(p1).rnvec;
                    else
                        parent(i).rnvec = population(p2).rnvec;
                    end
                else
                    parent(i).rnvec = population(p2).rnvec;
                end
            end
            % offspring generation
            count = count + 1;
            for i = 1:2:child_pop-1
                child(count) = Chromosome;
                child(count + 1) = Chromosome;
                p1 = i; p2 = i+1;
                [child(count).rnvec, child(count+1).rnvec] = ...
                    Evolve.crossover(parent(p1).rnvec, parent(p2).rnvec, muc, dim);
                child(count).rnvec = Evolve.mutate(child(count).rnvec, mum, dim);
                child(count+1).rnvec = Evolve.mutate(child(count+1).rnvec, mum, dim);
                count = count + 2;
            end
            if length(child) > pop
                child(pop+1:end) = [];
            end
        end
    else
        % parent pop selection
        for i = 1:pop
            parent(i) = Chromosome;
            p1 = 1+round(rand*(pop-1));
            p2 = 1+round(rand*(pop-1));
            if population(p1).front < population(p2).front
                parent(i).rnvec = population(p1).rnvec;
            elseif population(p1).front == population(p2).front
                if population(p1).CD > population(p2).CD
                    parent(i).rnvec = population(p1).rnvec;
                else
                    parent(i).rnvec = population(p2).rnvec;
                end
            else
                parent(i).rnvec = population(p2).rnvec;
            end
        end
        % offspring generation
        count = 1;
        for i = 1:2:pop-1
            child(count) = Chromosome;
            child(count + 1) = Chromosome;
            p1 = i; p2 = i+1;
            [child(count).rnvec, child(count+1).rnvec] = ...
                Evolve.crossover(parent(p1).rnvec, parent(p2).rnvec, muc, dim);
            child(count).rnvec = Evolve.mutate(child(count).rnvec, mum, dim);
            child(count+1).rnvec = Evolve.mutate(child(count+1).rnvec, mum, dim);
            count = count + 2;
        end
    end
    % evaluate offspring
    for i = 1:pop
        child(i) = evaluate(child(i), f, L, U);
    end
    clear intpopulation
    for i = 1:2*pop
        if i <= pop
            population(i) = reset(population(i));
            intpopulation(i) = population(i);
        else
            intpopulation(i) = child(i-pop);
        end
    end
    
    % non-dominant sorting
    [intpopulation, frontnumbers] = SolutionComparison.nondominatedsort(...
        intpopulation, 2*pop, no_of_objs);
    [intpopulation, minimums] = SolutionComparison.diversity(...
        intpopulation, frontnumbers, 2*pop, no_of_objs);
    population(1:pop) = intpopulation(1:pop);
    
    % convergence testing
    obj_data = vec2mat([population.objectives], no_of_objs);
    IGD = 0;
    if no_of_objs == 2
        for i = 1:size(data,1)
            c1 = data(i,1)*ones(pop,1);
            c2 = data(i,2)*ones(pop,1);
            IGD = IGD + sqrt(min(sum((obj_data-[c1 c2]).^2,2)));
        end
    elseif no_of_objs == 3
        for i = 1:size(data,1)
            c1 = data(i,1)*ones(pop,1);
            c2 = data(i,2)*ones(pop,1);
            c3 = data(i,3)*ones(pop,1);
            IGD = IGD + sqrt(min(sum((obj_data-[c1 c2 c3]).^2,2)));
        end
    else
        error('Too many objectives.')
    end
    store(generation)=IGD/size(data,1);  
    disp(['Generation ', num2str(generation), ' IGD score: ', num2str(store(generation))]);
end

end