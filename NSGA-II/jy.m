function mop=jy(testname,dimension)
%%run for zdt problems
mop=struct('name',[],'od',[],'pd',[],'domain',[],'func',[]);
switch lower(testname)
    case 'jy1'
        mop=jy1(mop,dimension);
    case 'jy2'
        mop=jy2(mop,dimension);
    case 'jy3'
        mop=jy3(mop,dimension);        
    case 'jy4'
        mop=jy4(mop,dimension);        
    case 'jy5'
        mop=jy5(mop,dimension);        
    case 'jy6'
        mop=jy6(mop,dimension);
    case 'jym4'
        mop=jym4(mop,dimension);        
    otherwise
        error('Undefined test problem name');
end
end

function p =jy1(p,pd)
 p.name     = 'JY1';
 p.pd      = pd;% default pd=30
 p.od   = 2;
% p.pf       = load('JY_F1.pf');
 p.domain   = [zeros(1,pd); ones(1,pd)];
 p.func     = @evaluate;
 
function F = evaluate(X)
    N   = pd;
    Y=sum((X(2:N)-sin(0.5*pi*X(2:N))).^2-cos(2*pi*(X(2:N)-sin(0.5*pi*X(2:N)))));
    G   = 2*sin(0.5*pi*X(1))*(N-1+Y);
    F   = zeros(1,2);

    F(1)= (1+G)*X(1);
    F(2)= (1+G)*(1-sqrt(X(1)))^5;
end
end

function p =jy2(p,pd)
 p.name     = 'JY2';
 p.pd       = pd;% default pd=30
 p.od       = 2;
% p.pf       = load('JY_F2.pf');
 p.domain   = [zeros(1,pd); ones(1,pd)];
 p.func     = @evaluate;
 
function F = evaluate(X)
    N   = pd;
    Y=sum((X(2:N)-sin(0.5*pi*X(2:N))).^2-cos(2*pi*(X(2:N)-sin(0.5*pi*X(2:N)))));
    G   = 2*sin(0.5*pi*X(1))*(N-1+Y);
    F   = zeros(1,2);

    F(1)= (1+G)*(1-X(1));
    F(2)= (1/2)*(1+G)*(X(1)+sqrt(X(1))*cos(4*pi*X(1))^2);
end
end

function p =jy3(p,pd)
 p.name     = 'JY3';
 p.pd       = pd;% default pd=30
 p.od       = 2;
%  p.pf       = load('JY_F3.pf');
 p.domain   = [zeros(1,pd); ones(1,pd)];
 p.func     = @evaluate;
 
function F = evaluate(X)
    N   = pd;
    Y=sum((X(2:N)-sin(0.5*pi*X(2:N))).^2-cos(2*pi*(X(2:N)-sin(0.5*pi*X(2:N)))));
    G   = 2*sin(0.5*pi*X(1))*(N-1+Y);
    F   = zeros(1,2);

    F(1)= (1+G)*X(1);
    F(2)= (1/2)*(1+G)*(1-X(1)^0.1+(1-sqrt(X(1)))^2*cos(3*pi*X(1))^2);
end
end

function p =jy4(p,pd)
 p.name     = 'JY4';
 p.pd       = pd;% default pd=30
 p.od       = 3;
%  p.pf       = load('JY_F4.pf');
 p.domain   = [ones(1,pd); 4*ones(1,pd)];
 p.func     = @evaluate;
 
function F = evaluate(X)
    N   = pd;
    G=sum((X(4:N)-2).^2);
    F   = zeros(1,3);

    F(1)= (1+G)*(X(1)/sqrt(X(2)*X(3)));
    F(2)= (1+G)*(X(2)/sqrt(X(1)*X(3)));
    F(3)= (1+G)*(X(3)/sqrt(X(1)*X(2)));
end
end

function p =jym4(p,pd)
 p.name     = 'JYm4';
 p.pd       = pd;% default pd=30
 p.od       = 3;
%  p.pf       = load('JY_mF4.pf');
 p.domain   = [ones(1,pd); 10*ones(1,pd)];
 p.func     = @evaluate;
 
function F = evaluate(X)
    N   = pd;
    G=sum((X(4:N)-5).^2);
    F   = zeros(1,3);

    F(1)= (1+G)*(X(1)/sqrt(X(2)*X(3)));
    F(2)= (1+G)*(X(2)/sqrt(X(1)*X(3)));
    F(3)= (1+G)*(X(3)/sqrt(X(1)*X(2)));
end
end

function p =jy5(p,pd)
 p.name     = 'JY5';
 p.pd       = pd;% default pd=30
 p.od       = 3;
%  p.pf       = load('JY_F5.pf');
 p.domain   = [zeros(1,pd); ones(1,pd)];
 p.func     = @evaluate;
 
function F = evaluate(X)
    N   = pd;
    G=sum((X(3:N)-0.5).^2);
    F   = zeros(1,3);

    F(1)= (1+G)*((1-X(1))*X(2));
    F(2)= (1+G)*((1-X(2))*X(1));
    F(3)= (1+G)*(1-X(1)-X(2)+2*X(1)*X(2))^6;
end
end

function p =jy6(p,pd)
 p.name     = 'JY6';
 p.pd       = pd;% default pd=30
 p.od       = 3;
%  p.pf       = load('JY_F6.pf');
 p.domain   = [zeros(1,pd); ones(1,pd)];
 p.func     = @evaluate;
 
function F = evaluate(X)
    N   = pd;
    G=(1/10)*(N-3+1+sum(X(3:N).^2-cos(2*pi*X(3:N))));
    F   = zeros(1,3);

    F(1)= (cos(0.5*pi*X(1))^4)*(cos(0.5*pi*X(2))^4);
    F(2)= (cos(0.5*pi*X(1))^4)*(sin(0.5*pi*X(2))^4);
    F(3)= ((1+G)/(1+cos(0.5*pi*X(1))^2))^(1/(1+G));
end
end