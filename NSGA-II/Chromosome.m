classdef Chromosome
    properties
        rnvec;
        objectives;
        convio;
        front;
        CD;
        dominationcount=0;
        dominatedset=[];
        dominatedsetlength=0;
    end
    methods
        function object = initialize(object,dim)
            object.rnvec=rand(1,dim);
        end
                
        function object = evaluate(object,f,L,U)
            var=L+object.rnvec.*(U-L);
%             dim=length(var);
            if ischar(f)
%                 if nargin == 5
%                     O = zeros(1,dim-1);
%                 end
%                 if strcmp(f,'ZDT4-S')
%                     var(2:dim) = 420.96 + M*(var(2:dim)-420.96)'; 
%                 elseif strcmp(f,'ZDT1')
%                     var(2:dim) = 1 + M*(var(2:dim)-1)';
%                 else
%                     var(2:dim) = M*(var(2:dim) - O)';
%                 end
                [object.objectives,object.convio] = mo_test_function(var,f);
            else
                object.convio = 0;
                object.objectives = evaluate(f, var');
                object.objectives = object.objectives';
            end
        end
        
        function object=reset(object)
            object.dominationcount=0;
            object.dominatedset=[];
            object.dominatedsetlength=0;
        end
    end
end