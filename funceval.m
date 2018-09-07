function fitness = funceval(population,dims,pop)
    % The trap-5 function
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

%     fitness = zeros(pop,1);
%     for i = 1:pop
%         fitness(i) = sum(population(i,1:10)) + 10 - sum(population(i,dims-9:dims));
%     end
end