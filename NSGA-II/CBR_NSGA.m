function store= CBR_NSGA(f,no_of_objs,L,U,pop)
load allmodels_CBR

injection = 0.1 * pop; % injection percentage
gen=100;
muc=20;% crossover index SBX
mum=10;% mutation index (polynomial mutation)
TrInt = 5; % transfer interval

if mod(pop,2)~=0
    pop=pop+1;
end
dim = length(L);

%     if nargin < 8
%         O = zeros(1,dim-1);
%         if nargin < 7
%             M = eye(dim-1);
%         end
%     end

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
    if mod(generation-2, TrInt)==0
        % inject solutions from case base
        % num_sol = size(allmodels_CBR, 1);
        bestRnvec = [population([population.front]==1).rnvec];
        bestRnvec = reshape(bestRnvec, [dim, length(bestRnvec)/dim])';
        % dist = sum((allmodels_CBR - bestRnvec).^2, 2);
        dist = sq_dist(allmodels_CBR(:,1:dim)', bestRnvec');
%                 [~, ind] = sort(dist);
        for i = 1:injection
%                     child(i).rnvec = allmodels_CBR(ind(i),:);
            [~,ind] = min(dist(:));
            [I_row, I_col] = ind2sub(size(dist),ind);
            dist(I_row,:) = inf;
            dist(:, I_col) = inf;
            child(i).rnvec = allmodels_CBR(I_row,1:dim);
        end
        count = injection + 1;
    else 
        count = 1;
    end
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
    for i=count:2:pop-1
        child(count)=Chromosome;
        child(count+1)=Chromosome;
        p1=i;
        p2=i+1;
        [child(count).rnvec,child(count+1).rnvec]=Evolve.crossover(parent(p1).rnvec,parent(p2).rnvec,muc,dim);
        child(count).rnvec = Evolve.mutate(child(count).rnvec,mum,dim);
        child(count+1).rnvec=Evolve.mutate(child(count+1).rnvec,mum,dim);
        count=count+2;
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
end  
        
%         figure(01)
%         hold on
%         plot(store)
%         
%         f1=[];
%         f2=[];
%         f3=[];      
%         population=population([population.front]==1);
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
%         figure()
%         if no_of_objs == 2
%             plot(f1,f2,'o')
%         else
%             plot3(f1,f2,f3,'o')
%         end    
        
end