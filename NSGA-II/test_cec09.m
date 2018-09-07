% test_cec09.m
% 
% Matlab source codes
% 
% Test function evaluation for CEC 2009 MOO Competition
%
% Usage: test_cec09(problem_name, vairaible_dimension), for example
% test_cec09('CF1',10)
% 
% Please refer to the report for more information.
%

function test_cec09(name, dim)

% get the search boundary
xrange = xboundary(name, dim);

NN = 10;
% randmoly generate NN points in the search space
x  = repmat(xrange(:,1),[1,NN]) + repmat((xrange(:,2)-xrange(:,1)),[1,NN]).*rand(dim, NN);

% %% === test matlab version with constraints
% % get the function handle
% fobj = cec09(name);
% % test the function
% [y,c] = fobj(x);
% % display the results
% disp('objectives');  disp(y);
% disp('constraints'); disp(c);

%% === test matlab version withot constraints
% get the function handle
fobj = cec09(name);
% test the function
y = fobj(x);
% display the results
disp('objectives');  disp(y);

% %% === test matlab-c version with constraints
% [y,c] = cec09m(x, name);
% disp('objectives');  disp(y);
% disp('constraints'); disp(c);

% %% === test matlab-c version without constraints
% y = cec09m(x, name);
% disp('objectives');  disp(y);

end