% generate_front_data.m
% 
% The Matlab source codes to generate the PF data of the test
%   instances for CEC 2009 Multiobjective Optimization Competition.
% 
% Usage: generate_front_data(), the result data files are in pf_data folder
% 
% Please refer to the report for more information.
% 
% History:
%   v1 Sept.08 2008

function generate_front_data()

PROBLEMS= ['UF1 '; 'UF2 '; 'UF3 '; 'UF4 '; 'UF5 '; 'UF6 '; 'UF7 '; 'UF8 '; 'UF9 '; 'UF10';];
DIMX    = [30 30 30 30 30 30 30 30 30 30];
% NOP     = [1000 1000 1000 1000 1000 1000 1000 10000 10000 10000];
NOP = 51*ones(1,10);
for p=1:10
    [PF,PS] = pareto( deblank(PROBLEMS(p,:)), NOP(p), DIMX(p) );
    PF      = PF';
    f       = sprintf('pf_data/%s.dat',deblank(PROBLEMS(p,:)));
    save(f, 'PF', '-ascii', '-tabs');
    clear PF PS;
end

PROBLEMS= ['CF1 '; 'CF2 '; 'CF3 '; 'CF4 '; 'CF5 '; 'CF6 '; 'CF7 '; 'CF8 '; 'CF9 '; 'CF10';];
DIMX    = [10 10 10 10 10 10 10 10 10 10];
NOP     = [1000 1000 1000 1000 1000 1000 1000 10000 10000 10000]; 
for p=1:10
    [PF,PS] = pareto( deblank(PROBLEMS(p,:)), NOP(p), DIMX(p) );
    PF      = PF';
    f       = sprintf('pf_data/%s.dat',deblank(PROBLEMS(p,:)));
    save(f, 'PF', '-ascii', '-tabs');
    clear PF PS;
end

end