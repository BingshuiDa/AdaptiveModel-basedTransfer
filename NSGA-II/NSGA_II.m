function store= NSGA_II(f,no_of_objs,L,U,pop)
    
load allmodels_TS
load allmodels_CBR
load allmodels_AE
global maxdim;

gen=250;
muc=15;% crossover index SBX
mum=20;% mutation index (polynomial mutation)

if mod(pop,2)~=0
    pop=pop+1;
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

solutions_history = zeros(pop, dim, gen);
    
for i =1:pop
    population(i)=Chromosome;
    population(i)=initialize(population(i),dim);
end
for i=1:pop
	population(i)=evaluate(population(i),f,L,U);
end 
[population,frontnumbers]=SolutionComparison.nondominatedsort(population,pop,no_of_objs);
[population,minimums]=SolutionComparison.diversity(population,frontnumbers,pop,no_of_objs);
%         [population,minimums]=SolutionComparison.extremumbiasing(population,frontnumbers,pop,no_of_objs);
for generation=1:gen
    for i = 1:pop
        parent(i)=Chromosome();
        p1=1+round(rand(1)*(pop-1));
        p2=1+round(rand(1)*(pop-1));
        if population(p1).front < population(p2).front
            parent(i).rnvec=population(p1).rnvec;
        elseif population(p1).front == population(p2).front
            if population(p1).CD > population(p2).CD
                parent(i).rnvec=population(p1).rnvec;
            else 
                parent(i).rnvec=population(p2).rnvec;
            end
        else
        	parent(i).rnvec=population(p2).rnvec;
        end
    end
    solutions = zeros(pop, dim);
    for i = 1:pop
        solutions(i,:) = parent(i).rnvec;
    end
    count=1;
    for i=1:2:pop-1
        child(count)=Chromosome;
        child(count+1)=Chromosome;
        p1=i;
        p2=i+1;
        [child(count).rnvec,child(count+1).rnvec]=Evolve.crossover(parent(p1).rnvec,parent(p2).rnvec,muc,dim);
        child(count).rnvec = Evolve.mutate(child(count).rnvec,mum,dim);
        child(count+1).rnvec=Evolve.mutate(child(count+1).rnvec,mum,dim);
        count=count+2;
    end

    for i = 1:pop
        child(i).rnvec(child(i).rnvec < 0) = 0;
        child(i).rnvec(child(i).rnvec > 1) = 1;
    end

    for i=1:pop
        child(i)=evaluate(child(i),f,L,U);
    end
    clear intpopulation
    for i=1:2*pop
        if i <= pop
            population(i)=reset(population(i));
            intpopulation(i)=population(i);
        else
            intpopulation(i)=child(i-pop);
        end
    end

    [intpopulation,frontnumbers]=SolutionComparison.nondominatedsort(intpopulation,2*pop,no_of_objs);
    [intpopulation,minimums]=SolutionComparison.diversity(intpopulation,frontnumbers,2*pop,no_of_objs);           
    %            [intpopulation,minimums]=SolutionComparison.extremumbiasing(intpopulation,frontnumbers,2*pop,no_of_objs);
    population(1:pop)=intpopulation(1:pop);

%            Convergence testing
    obj_data = vec2mat([population.objectives],no_of_objs);           
    IGD = 0;
    if mod(generation,100) == 0
    	disp(generation)
    end
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
        error('Too many objectives.');
    end
    store(generation)=IGD/size(data,1);  
    disp(['Generation ', num2str(generation), ' IGD score: ', num2str(store(generation))]);

%     solutions = [];
    solutions = [];
    for i = 1:length(population)
%         if population(i).front == 1
%             solutions = [solutions; population(i).rnvec];
%         end
        solutions = [solutions; population(i).rnvec];
    end
    solutions = solutions + 0.05*randn(size(solutions)); % corrupted by noise
    solutions(solutions > 1) = 1;
    solutions(solutions < 0) = 0;
    if dim ~= maxdim
        solutions = [solutions, rand(size(solutions, 1), maxdim - dim)];
    end
    allmodels_CBR = [allmodels_CBR; solutions];
    solutions_history(:, :, generation) = solutions;

%            f1=[];
%            f2=[];
%            for i=1:length(population)
%               f1=[f1,population(i).objectives(1)];
%               f2=[f2,population(i).objectives(2)];
%            end   
%            plot(f1,f2,'o')
%            drawnow
end  

% save data for autoencoder
sourceAE.history = solutions_history;
solutions = solutions(:, 1:dim);
sourceAE.best = solutions;
allmodels_AE{length(allmodels_AE)+1} = sourceAE;
save('allmodels_AE', 'allmodels_AE');

% build probabilistic model for Transfer Stacking
% pop_model = round(1.0*pop);
pop_model = pop;
solutions = zeros(pop_model, dim);
for i = 1:pop_model
    if i <= pop
        solutions(i,:) = population(i).rnvec;
    else
        solutions(i,:) = rand(1,dim);
    end
end
if dim ~= maxdim
	solutions = [solutions, rand(size(solutions, 1), maxdim - dim)];
end
model = ProbabilityModel('mvarnorm');
model = ProbabilityModel.buildmodel(model, solutions);
allmodels_TS{length(allmodels_TS)+1} = model;
save('allmodels_TS','allmodels_TS');

save('allmodels_CBR', 'allmodels_CBR');


%         f1=[];
%         f2=[];
%         f3=[];      
%         population=population([population.front]==1);
        
        % store pareto optimal solutions in case base for CBR
%         solutions = zeros(length(population), dim);
%         for i = 1:length(population)
%             solutions(i,:) = population(i).rnvec;
%         end
%         allmodels_CBR = [allmodels_CBR; solutions];
        
%         for i=1:length(population)
%             if no_of_objs == 2
%                 f1=[f1,population(i).objectives(1)];
%                 f2=[f2,population(i).objectives(2)];
%             else
%                 f1=[f1,population(i).objectives(1)];
%                 f2=[f2,population(i).objectives(2)];
%                 f3=[f3,population(i).objectives(3)];
%             end
%         end   
%         if no_of_objs == 2
%             plot(f1,f2,'o')
%         else
%             plot3(f1,f2,f3,'o')
%         end    
%         figure(02)
%         plot(store)
end