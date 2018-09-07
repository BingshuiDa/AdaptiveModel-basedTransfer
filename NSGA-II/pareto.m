% pareto.m
% 
% The Matlab source codes to generate the PF and the PS of the test
%   instances for CEC 2009 Multiobjective Optimization Competition.
% 
% Usage: [pf, ps] = pareto(problem_name, no_of_points, variable_dim)
% 
% Please refer to the report for more information.
% 
% History:
%   v1 Sept.08 2008

function [pf, ps] = pareto(name, no, dim)

    if nargin<3, dim = 3; end
    if nargin<2, no  = 500; end
    switch name
        case 'UF1'
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-sqrt(pf(1,:));
            ps          = zeros(dim,no);
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
        case 'UF2'
            pf            = zeros(2,no);
            pf(1,:)       = linspace(0,1,no);
            pf(2,:)       = 1-sqrt(pf(1,:));
            ps            = zeros(dim,no);
            ps(1,:)       = linspace(0,1,no);
            X1            = repmat(ps(1,:),[dim,1]);
            A             = 6.0*pi*X1 + pi/dim*repmat((1:dim)',[1,no]);
            B             = 0.3*X1.*(X1.*cos(4.0*A)+2.0).*cos(A);
            ps(3:2:dim,:) = B(3:2:dim,:);
            B             = 0.3*X1.*(X1.*cos(4.0*A)+2.0).*sin(A);
            ps(2:2:dim,:) = B(2:2:dim,:);
            clear A B;
        case 'UF3'
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-sqrt(pf(1,:));
            ps          = zeros(dim,no);
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = repmat(ps(1,:),[dim-1,1]).^(0.5+1.5*(repmat((0:1:(dim-2))',[1,no]))/(dim-2.0));
        case 'UF4'
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-pf(1,:).^2;
            ps          = zeros(dim,no);
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
        case 'UF5'
            no          = 21;
            pf          = zeros(2,no);
            pf(1,:)     = (0:1:20)/20.0;
            pf(2,:)     = 1-pf(1,:);
            ps          = zeros(dim,no);
            ps(1,:)     = pf(1,:);
            ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]));
        case 'UF6'
            num                     = floor(no/3);
            pf                      = zeros(2,no);
            pf(1,1:num)             = 0.0;
            pf(1,(num+1):(2*num))   = linspace(0.25,0.5,num);
            pf(1,(2*num+1):no)      = linspace(0.75,1.0,no-2*num);
            pf(2,:)                 = 1-pf(1,:);
            ps                      = zeros(dim,no);
            ps(1,:)                 = pf(1,:);
            ps(2:dim,:)             = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
        case 'UF7'
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-pf(1,:);
            ps          = zeros(dim,no);
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
        case {'UF8','UF10'}
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            pf          = zeros(3,no);
            pf(1,:)     = cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:));
            pf(2,:)     = cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:));
            pf(3,:)     = sin(0.5*pi*ps(1,:));   
            clear s t;
        case 'UF9'
            num             = floor(sqrt(no));
            no              = num*num;
            noA             = floor(num/2);
            A               = zeros(1,num);
            A(1,1:noA)      = linspace(0,0.25,noA);
            A(1,noA+1:num)  = linspace(0.75,1,num-noA);
            [s,t]           = meshgrid(A,linspace(0,1,num));
            ps              = zeros(dim,no);
            ps(1,:)         = reshape(s,[1,no]);
            ps(2,:)         = reshape(t,[1,no]);            
            ps(3:dim,:)     = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            pf              = zeros(3,no);
            pf(1,:)         = ps(1,:).*ps(2,:);
            pf(2,:)         = (1.0-ps(1,:)).*ps(2,:);
            pf(3,:)         = 1.0-ps(2,:);    
            clear A s t;
        case 'CF1'
            no          = 21;
            pf          = zeros(2,no);
            pf(1,:)     = (0:1:20)/20.0;
            pf(2,:)     = 1-pf(1,:);
            ps          = [];
        case 'CF2'
            no          = floor(no/3.0);
            pf          = zeros(2,3*no);
            pf(1,(no+1):(2*no))   = linspace(0.25^2,0.25,no);
            pf(1,(2*no+1):(3*no)) = linspace(0.75^2,1.0,no);
            pf(2,:)     = 1-sqrt(pf(1,:));
            ps          = []; 
        case 'CF3'
            no          = floor(no/3.0);
            pf          = zeros(2,3*no);
            pf(1,(no+1):(2*no))     = linspace(0.25^0.5,0.5^0.5,no);
            pf(1,(2*no+1):(3*no))   = linspace(0.75^0.5,1.0^0.5,no);
            pf(2,:)     = 1-pf(1,:).^2;
            ps          = [];
        case {'CF4','CF5'}
            no          = floor(no/4.0);
            pf          = zeros(2,4*no);
            pf(1,1:(2*no))          = linspace(0,0.5,2*no);
            pf(2,1:(2*no))          = 1-pf(1,1:(2*no));
            pf(1,(2*no+1):(3*no))   = linspace(0.5,0.75,no);
            pf(2,(2*no+1):(3*no))   = -0.5*pf(1,(2*no+1):(3*no))+0.75;
            pf(1,(3*no+1):(4*no))   = linspace(0.75,1.0,no);
            pf(2,(3*no+1):(4*no))   = 1.125-pf(1,(3*no+1):(4*no));
            ps          = [];
        case {'CF6','CF7'}
            no          = floor(no/4.0);
            pf          = zeros(2,4*no);
            pf(1,1:(2*no))          = linspace(0,0.5,2*no);
            pf(2,1:(2*no))          = (1.0-pf(1,1:(2*no))).^2;
            pf(1,(2*no+1):(3*no))   = linspace(0.5,0.75,no);
            pf(2,(2*no+1):(3*no))   = 0.5*(1.0-pf(1,(2*no+1):(3*no)));
            pf(1,(3*no+1):(4*no))   = linspace(0.75,1.0,no);
            pf(2,(3*no+1):(4*no))   = 0.25*sqrt(1.0-pf(1,(3*no+1):(4*no)));
            ps          = []; 
        case 'CF8'
            no          = floor(no/5.0);
            pf          = zeros(3,5*no);
            for k = 0:1:4
                s = k*no+1; 
                e = k*no+no;
                pf(3,s:e) = linspace(0,1,no);
                pf(1,s:e) = sqrt(k/4.0*(1.0-pf(3,s:e).^2));
                pf(2,s:e) = sqrt(1-pf(1,s:e).^2-pf(3,s:e).^2);
            end
            ps          = []; 
        case {'CF9','CF10'}
            no1         = floor(sqrt(no/12));
            no          = 4*no1*no1;
            pf          = zeros(3,3*no);
            s           = 1; 
            e           = no;
            pf(1,s:e)   = 0;
            pf(2,s:e)   = linspace(0,1,no);
            pf(3,s:e)   = sqrt(1-pf(2,s:e).^2);            
            for k = 1:1:2
                s = k*no+1; 
                e = k*no+no;
                A = repmat(linspace(0,1,4*no1),[no1,1]);
                B = repmat((linspace(2*k-1, 2*k, no1)/4.0)',[1,4*no1]);
                B = sqrt(B.*repmat(1-linspace(0,1,4*no1).^2,[no1,1]));
                pf(3,s:e) = reshape(A,1,no);
                pf(1,s:e) = reshape(B,1,no);
                pf(2,s:e) = sqrt((1-pf(1,s:e).^2-pf(3,s:e).^2));
            end
            ps          = [];             
    end
end
