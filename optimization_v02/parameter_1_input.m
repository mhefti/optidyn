clc; close all; mout = 'matlaboutput.txt';
% clear all
if exist(mout,'file')
    fid = fopen(mout,'w+');
    fclose(fid);
end
diary(mout);
%% input
exp_list = [140714 140718 140730]; % choose dates of expts
timecut_list = 3600*[7 9 5]; % choose at what time to cut the expt
mode_list = ['a' 'a' 'a']; % choose also the mode (a: ads/d: des) 

% choose what to fit
fit_y = true; % composition of H2O: yH2O
fit_T = true; % temperature at position 0.5*L

% choose to fit on brutus or not
brutmode = 'nobrutus'; % or brutus/nobrutus

% if brutus, choose number of workers (local machine: default 2 workers)
slaves = 16;

% choose the algorithm to be used
algorithm = 'GlobalSearch'; 
parallel = 'noparallel'; % or 'noparallel'
% available options: 
% -----------------------------
% algorithm         parallel?       comments:
% lsq               no
% GlobalSearch      no              fmincon/sqp
% MultiStart        yes             fmincon/sqp
% patternsearch     yes
% -----------------------------

if strcmp(brutmode,'brutus')
    brutus = true;
else
    brutus = false;
end

%% set up experiments output
numexp = 1:length(exp_list);
expinfo.list = exp_list;
expinfo.numexp = max(numexp);
expinfo.modes = mode_list;

condmatall = [];
% structure-fieldname function, % e.g. data.exp_140714 = ...
% valid inputs ('yH2O',i), ('T05',i)
sfield = @(q,ii) sprintf('%s_%s_%s',q,num2str(expinfo.list(ii)),mode_list(ii));

% save the relevant data in data.('expdate')
for i = numexp
    data.(sfield('exp',i)) = get_experiment_data(exp_list(i),...
        brutmode,timecut_list(i),mode_list(i));
    condmatall = [condmatall data.(sfield('exp',i)).condmat];
end
close all;
% write conditions
dlmwrite('prms/conditions.dat',condmatall,'\t')
% write iterations - 
fid = fopen('prms/iterations.dat','w+');
tail = '// first/last experiment to evaluate';
fprintf(fid,'%s',sprintf('%i %i\t\t%s',1,expinfo.numexp,tail));
fclose(fid);

expinfo.conditions = condmatall;
% expinfo.fparamids = [16 18]; % rows of mtc and htc 
expinfo.fparamids = [16 12];
expinfo.fity = fit_y;
expinfo.fitT = fit_T;

%% set up fitting
fid = fopen('parameters/initial_guess.txt','r');
initvals = textscan(fid,'%f','delimiter','','HeaderLines',0);
initvals = initvals{1}(:);
fclose(fid);

% experimental data
for i = numexp
    
    if fit_y && fit_T
        y.(sfield('yH2O',i)) = data.(sfield('exp',i)).yH2O;
        y.(sfield('T05',i)) = data.(sfield('exp',i)).T05 + 273.15;
        expinfo.numpoints(i) = length(data.(sfield('exp',i)).yH2O);
        expinfo.maxy(i) = max(y.(sfield('yH2O',i)));
        expinfo.maxT(i) = max(y.(sfield('T05',i)));
    elseif fit_y
        y.(sfield('yH2O',i)) = data.(sfield('exp',i)).yH2O;
        expinfo.numpoints(i) = length(data.(sfield('exp',i)).yH2O);
        expinfo.maxy(i) = max(y.(sfield('yH2O',i)));
    else
        y.(sfield('T05',i)) = data.(sfield('exp',i)).T05 + 273.15;
        expinfo.numpoints(i) = length(data.(sfield('exp',i)).T05);
        expinfo.maxT(i) = max(y.(sfield('T05',i)));
    end
    
end
x = [];

% constraining
lb = 1e0*ones(size(initvals));
lb(1) = 1e-3;
lb(2) = 1e-3;
lb(4) = 1e-3;
lb(3) = 1;
lb(5) = 2;
lb(6) = 0.001;
lb(7) = 10000;

ub = 1e4*ones(size(initvals));
ub(1) = 5e3;
ub(2) = 5e3;
ub(3) = 25;
ub(4) = 500;
ub(5) = 20;
ub(6) = 0.1;
ub(7) = 45000;

% exectuable name
execname = 'hc_v02_J'; % use short names for brutus

% initialize log-file
numpars = length(initvals);
param_number_func = @(ii) sprintf('parameter_%2.2d\t',ii);
headerfunc = @(ii) sprintf('%s resnorm\t',param_number_func(ii));

flog = 'fitting_log.txt';
fid = fopen(flog,'w');
fprintf(fid,'%s\n',headerfunc(1:numpars));
fclose(fid);

if strcmp(algorithm,'patternsearch') ~= true
    ams = {'interior-point','trust-region-reflective','sqp','active-set'};
    options = optimset('Display','final','Algorithm', ams{3},...
    'TolX',1e-4,'TolFun',1e-4, 'DiffMinChange', 0.01,'PlotFcns',...
    @optimplotfval);
else % options for patternsearch
    options = psoptimset();
    options.Cache = 'on';
    options.TolMesh = 1e-6;
    options.TolX = 1e-4;
    options.TolFun = 1e-4;
	options.MaxIter = 1e3*length(initvals);
	options.MaxFunEvals = 4000*length(initvals);
    options.Display = 'final';
    options.MeshAccelerator = 'on';
    optsions.PlotFcns = @psplotbestf;
end

tic;
% call the optimizer
optimoutput = getfit(@fitness_parallel,initvals,x,y,lb,ub,options,execname,...
        algorithm,parallel,brutmode,expinfo,slaves);
toc;

% write results
headerfunc = @(ii) sprintf('%s resnorm\t outputflag',param_number_func(ii));
numcols = numpars + 2;
fid = fopen('output.txt','w');
outrowvec = [optimoutput.params'.*initvals...
    optimoutput.resnorm optimoutput.flag];
fprintf(fid,'%s\n',headerfunc(1:numpars));
fprintf(fid,repmat('%f\t',1,numcols),outrowvec);
fclose(fid);