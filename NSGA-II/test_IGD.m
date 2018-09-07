% test_IGD.m
% 
% Matlab source codes
% 
% Test IGD metric for CEC 2009 MOO Competition
%
% Usage: test_IGD()
% 
% Please refer to the report for more information.
%

function test_IGD()

% use mex -setup to choose a compiler if necessary
disp('----------------------------------------');
disp('build IGD metric');
% mex IGD.cpp;
disp('----------------------------------------');

disp('----------------------------------------');
disp('test IGD without constraint on UF1')
% load the PF* data
PFStar  = load('pf_data/UF1.dat');
PFStar  = PFStar';

% randomly generate 1000 points inside the search space
xrange = xboundary('UF1', 30);
x  = repmat(xrange(:,1),[1,100]) + repmat((xrange(:,2)-xrange(:,1)),[1,100]).*rand(30, 100);
% evaluate the points
fobj = cec09('UF1');
PF   = fobj(x);


str  = sprintf('%.5f = IGD(PFStar, PF)', IGD(PFStar, PF));
disp(str);
clear PFStar xrange x fobj PF;
disp('----------------------------------------');

disp('----------------------------------------');
disp('test IGD with constraint on CF1')
% load the PF* data
PFStar  = load('pf_data/CF1.dat');
PFStar  = PFStar';

% randomly generate 1000 points inside the search space
xrange = xboundary('CF1', 10);
x      = repmat(xrange(:,1),[1,100]) + repmat((xrange(:,2)-xrange(:,1)),[1,100]).*rand(10, 100);
% evaluate the points
fobj   = cec09('CF1');
[PF,C] = fobj(x);
str  = sprintf('%.5f = IGD(PFStar, PF, C)', IGD(PFStar, PF, C));
disp(str);
clear PFStar xrange x fobj PF C;
disp('----------------------------------------');

end