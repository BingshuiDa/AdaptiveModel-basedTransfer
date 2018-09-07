function mop=glt(testname,dimension)
%%run for zdt problems
mop=struct('name',[],'od',[],'pd',[],'domain',[],'func',[]);
switch lower(testname)
    case 'glt1'
        mop=glt1(mop,dimension);
    case 'glt2'
        mop=glt2(mop,dimension);
    case 'glt3'
        mop=glt3(mop,dimension);        
    case 'glt4'
        mop=glt4(mop,dimension);        
    case 'glt5'
        mop=glt5(mop,dimension);        
    case 'glt6'
        mop=glt6(mop,dimension);                    
    otherwise
        error('Undefined test problem name');
end
end
%%
function prob=glt1(prob,dim)
prob.name='GLT1';
prob.pd=dim;% default dim=10
prob.od=2;
prob.domain=[0,1;-1*ones(dim-1,1),ones(dim-1,1)];
% prob.pf=load('SMEAGLT1.pf');
prob.func=@evaluate;
    function y = evaluate(x)
        n         = size(x,1);
        g         = x-sin(2*pi*x(1)+(1:1:n)'*pi/n);
        y         = zeros(2,1);
        y(1)      = (1+sum(g).^2.0)*x(1);
        y(2)      = (1+sum(g).^2.0)*(2-x(1)-sign(cos(2*pi*x(1))));
    end
end
%%
function prob=glt2(prob,dim)
prob.name='GLT2';
prob.pd=dim;% default dim=10
prob.od=2;
prob.domain=[0,1;-1*ones(dim-1,1),ones(dim-1,1)];
% prob.pf=load('SMEAGLT2.pf');
prob.func=@evaluate;

    function y = evaluate(x)
        n         = size(x,1);
        g         = x-sin(2*pi*x(1)+(1:1:n)'*pi/n);
        y         = zeros(2,1);
        y(1)      = (1+sum(g).^2.0)*(1-cos(pi*x(1)/2));
        y(2)      = (1+sum(g).^2.0)*(10-10*sin(pi*x(1)/2));
    end
end
%%
function prob=glt3(prob,dim)
prob.name='GLT3';
prob.pd=dim;% default dim=10
prob.od=2;
prob.domain=[0,1;-1*ones(dim-1,1),ones(dim-1,1)];
% prob.pf=load('SMEAGLT3.pf');
prob.func=@evaluate;

    function y = evaluate(x)
        n         = size(x,1);
        g         = x-sin(2*pi*x(1)+(1:1:n)'*pi/n);
        y         = zeros(2,1);
        y(1)      = (1+sum(g).^2.0)*x(1);
        if y(1)<=0.05
            y(2)      = (1+sum(g).^2.0)*(1-19*x(1));
        else
            y(2)      = (1+sum(g).^2.0)*(1/19-x(1)/19);
        end
    end
end
%%
function prob=glt4(prob,dim)
prob.name='GLT4';
prob.pd=dim;% default dim=10
prob.od=2;
prob.domain=[0,1;-1*ones(dim-1,1),ones(dim-1,1)];
% prob.pf=load('SMEAGLT4.pf');
prob.func=@evaluate;

    function y = evaluate(x)
        n         = size(x,1);
        g         = x-sin(2*pi*x(1)+(1:1:n)'*pi/n);
        y         = zeros(2,1);
        y(1)      = (1+sum(g).^2.0)*x(1);
        y(2)      = (1+sum(g).^2.0)*(2-2*x(1).^(0.5)*cos(2*x(1).^0.5*pi).^2);
    end
end
%%
function prob=glt5(prob,dim)
prob.name='GLT5';
prob.pd=dim;% default dim=10
prob.od=3;
prob.domain=[0,1;0,1;-1*ones(dim-2,1),ones(dim-2,1)];
% prob.pf=load('SMEAGLT5.pf');
prob.func=@evaluate;

    function y = evaluate(x)
        n         = size(x,1);
        g         = x-sin(2*pi*x(1)+(1:1:n)'*pi/n);
        y         = zeros(3,1);
        y(1)      = (1+sum(g).^2.0)*(1-cos(x(1)*pi/2))*(1-cos(x(2)*pi/2));
        y(2)      = (1+sum(g).^2.0)*(1-cos(x(1)*pi/2))*(1-sin(x(2)*pi/2));
        y(3)      = (1+sum(g).^2.0)*(1-sin(x(1)*pi/2));
    end
end
%%
function prob=glt6(prob,dim)
prob.name='GLT6';
prob.pd=dim;% default dim=10
prob.od=3;
prob.domain=[0,1;0,1;-1*ones(dim-2,1),ones(dim-2,1)];
% prob.pf=load('SMEAGLT6.pf');
prob.func=@evaluate;

    function y = evaluate(x)
        n         = size(x,1);
        g         = x-sin(2*pi*x(1)+(1:1:n)'*pi/n);
        y         = zeros(3,1);
        y(1)      = (1+sum(g).^2.0)*(1-cos(x(1)*pi/2))*(1-cos(x(2)*pi/2));
        y(2)      = (1+sum(g).^2.0)*(1-cos(x(1)*pi/2))*(1-sin(x(2)*pi/2));
        y(3)      = (1+sum(g).^2.0)*(2-sin(x(1)*pi/2)-sign(cos(4*pi*x(1))));
    end
end