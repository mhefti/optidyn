close all; mout = 'matlaboutput.txt'; %clc;
% clear all
if exist(mout,'file')
    fid = fopen(mout,'w+');
    fclose(fid);
end
diary(mout);
%% input

% choose experiments to be fitted 
% -------------------------------------------------------------------------
% all experiments
% exp_list = [140801 140801 140808 140808 140812 140714 140718]; % choose dates of expts
% timecut_list = 3600*[0.8 1.5 0.6 1 5 7 9]; % choose at what time to cut the expt
% mode_list = ['a' 'd' 'a' 'd' 'a' 'a' 'a']; % choose also the mode (a: ads/d: des) 

exp_list = [140714 140718 140730]; % choose dates of expts
timecut_list = 3600*[7 9 5]; % choose at what time to cut the expt in [s]
mode_list = ['a' 'a' 'a']; % choose also the mode (a: ads/d: des) 




% objective variables 
% -------------------------------------------------------------------------
fit_y = true; % composition of H2O: yH2O
fit_T = true; % temperature at position 0.5*L

% fitting parameters
% -------------------------------------------------------------------------

% isotherm related parameters ---------------------------------------------
fitisotype = 'Sips_Sips';           % isotherm, 'Sips_Sips or 'Do_modified'

fit_iso = false;                    % if true, dynamic isotherm fitting to
                                    % static data then the scaling factor

numisopar = 6;                      % # of isotherm parameters
                                    % ** NOTE: must also be set to the
                                    % right value when fit_iso is enabled

isoflex = true;                     % flag for individual iso-par fitting
fitisolocs = [4:6];                 % provide which parameters of the 
                                    % isotherm should be fitted; don't have
                                    % to be subsequent: [1,3,5] possible
                                    % default: all, i.e. 1:numisopar
                                    % *** NOTE: can be combined with the
                                    % feature fit_iso 
                                    % *** CAUTION: if isoflex == true make
                                    % sure to have the right isotherm 
                                    % parameter values in prms/isotherm.dat 
                                    
% mass transfer parameters ------------------------------------------------
mtc_0 = 0.0035;                     % enter the base mass transfer coeff
                                    % *** NOTE: is overwritten if mtc0 is
                                    % fitted as well
mtcmodel = 'monotonic';             % options are: 'constant', 'tangens', 
                                    % 'monotonic', 'exponential'
fitmtc = false;                     % mass transfer coefficient
fitmtcmodel = true;                 % fit a mtc model, note that this
                                    % doesn't fit the mtc itself, but  
                                    % parameters of the scaling function
                                    % specified in mtcmodel
nummtcpar = 2;                      % number of mtcmodel parameters
                                    
% heat related parameters -------------------------------------------------
fitisothermal = false;              % hads will be = 0 in parameter1.dat
                                    % htc will be 1e6 in conditions.dat
fithtc = false;                     % heat transfer coefficient
fitHads = true;                     % heat of adsorption; note that Hads 
                                    % should be entered in (- kJ/mol)
% heat of adsorption scaling
hscaling = false;                    % flag if hads is scaled
hmodel = 'custom';                  % options are: 'custom'; any other
                                    % string results in the use of a
                                    % constant heat of adsorption for the
                                    % water component
numhscalepar = 3;                   % number of parameters for the scaling 

% -------------------------------------------------------------------------
% keep the order of the fitting vector, i.e. arrange as follows: 
% [     isopar(1)                   % 1
%          .                        % .
%          .                        % .
%       isopar(numisopar)           % numisopar
%       mtc                         % numisopar + 1
%       mtcmodelpar_1               % numisopar + 2
%       .                           % .
%       mtcmodelpar_X               % numisopar + nummtcpar + 1
%       htc                         % numisopar + nummtcpar + 2
%       Hads                        % numisopar + nummtcpar + 3
%       hscalepar_1                 % numisopar + nummtcpar + 4
%       .                           % .
%       hscalepar_Y              ]  % numpar
% -------------------------------------------------------------------------

% choose to fit on local machin, brutus or euler
brutmode = 'nobrutus'; %  'brutus' or 'nobrutus'
eulermode = true;      % *** NOTE: to run on euler, brutmode = 'brutus' !

% objective function type
phitype = 'DV';              % options: - maximum likelihood estimate 'MLE'
                             %          - derivative weighted: 'DV'

% if brutus, choose number of workers (local machine: default 2 workers)
slaves = 16;

% set the time for the jobs [h]
timespan = 36;

% set the timeout [s]
timeout = 300;

% choose the algorithm to be used
algorithm = 'GlobalSearch'; 
parallel = 'noparallel'; % or 'noparallel'

% optional for MultiStart, set number of initial points (default: 2e4)
numMSpoints = 2e4;

% available options: 
% -------------------------------------------------------------------------
% algorithm         parallel?       comments:
% lsq               no
% GlobalSearch      no              fmincon/sqp
% fminsearchbnd     no              derivative free
% MultiStart        yes             fmincon/sqp
% patternsearch     yes
% -------------------------------------------------------------------------

if strcmp(brutmode,'brutus')
    brutus = true;
else
    brutus = false;
end

% end of input
%% set up experiments output
numexp = 1:length(exp_list);
expinfo.list = exp_list;
expinfo.numexp = max(numexp);
expinfo.modes = mode_list;
expinfo.time = timespan;
expinfo.fitiso = fit_iso;
expinfo.time_out = timeout;
expinfo.fit_mtc = fitmtc;
expinfo.fit_htc = fithtc;
expinfo.fit_Hads = fitHads;
expinfo.fit_isothermal = fitisothermal;
expinfo.num_isopar = numisopar;
expinfo.iso_type = fitisotype;
expinfo.isoflex = isoflex;
expinfo.fitisolocs = fitisolocs;
expinfo.numMSpoints = numMSpoints; 
expinfo.fitmtcmodel = fitmtcmodel;
expinfo.mtcmodel = mtcmodel;
expinfo.mtcmodelid = 16;
expinfo.nummtcpar = nummtcpar;
expinfo.phitype = phitype;
expinfo.eulermode = eulermode;
expinfo.hscaling = hscaling; 
expinfo.numhscalepar = numhscalepar;
expinfo.hmodel = hmodel;
expinfo.hmodelid = 18;

% extract # of fitting parameters: 
numpar = numisopar;

if expinfo.fitiso
    numpar = numpar - numisopar + 1;
end

if ~expinfo.fitiso && expinfo.isoflex
    numpar = numel(expinfo.fitisolocs);
end

if expinfo.fit_mtc
    numpar = numpar + 1;
end

if expinfo.fit_htc
    numpar = numpar + 1;
end

if expinfo.fit_Hads
    numpar = numpar + 1;
end

if expinfo.fit_isothermal && expinfo.fit_Hads
    error('error: are you sure to fit Hads to isothermal conditions?')
end

if expinfo.fit_htc && expinfo.fit_htc 
    error('error: htc is to be fitted in isothermal conditions')
end

if expinfo.fitmtcmodel 
    numpar = numpar + nummtcpar;
end

if expinfo.hscaling
    numpar = numpar + numhscalepar;
end

expinfo.num_pars = numpar;

fprintf('the number of fitting parameters is:           %i\n',expinfo.num_pars)

condmatall = [];

% sfield constructs a string 'q_DATE_mode' where DATE is the date of the
% experiment specified in expinfo.list and mode is adsorption (a) or
% desorption (d) specified in expinfo.modes. 
% structure-fieldname function, % e.g. data.exp_140714_a = ...
% to be used as a field name; if so enclose by round brackets
% valid inputs ('yH2O',i), ('T05',i)
sfield = @(q,ii) sprintf('%s_%s_%s',q,num2str(expinfo.list(ii)),...
    mode_list(ii));

% save the relevant data in data.('expdate')
for i = numexp
    data.(sfield('exp',i)) = get_experiment_data(exp_list(i),...
        brutmode,timecut_list(i),mode_list(i),mtc_0);
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
expinfo.htc_id = 18;
expinfo.mtc_id = 16;
expinfo.Hads_id = 12;
expinfo.fity = fit_y;
expinfo.fitT = fit_T;

%% set up fitting
fid = fopen('parameters/initial_guess.txt','r');
ins = textscan(fid,'%f%f%f','delimiter','\t','HeaderLines',1);

lb = ins{1}(:);
initvals = ins{2}(:);
ub = ins{3}(:);

fclose(fid);

% experimental data
for i = numexp
    
    % save the time no matter what
    y.(sfield('time',i)) = data.(sfield('exp',i)).t;    
    
    % adapt the fields depending on the fitting flags
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
    else % only fit_T is true
        y.(sfield('T05',i)) = data.(sfield('exp',i)).T05 + 273.15;
        expinfo.numpoints(i) = length(data.(sfield('exp',i)).T05);
        expinfo.maxT(i) = max(y.(sfield('T05',i)));
    end
    
end

x = [];
% constraining
if length(initvals) ~= expinfo.num_pars
    fprintf('initial guesses:                               %i parameters\n',length(initvals))
    fprintf('to be fitted:                                  %i parameters\n',expinfo.num_pars)
    error('initial guess does not match the number of fitting parameters')
end

% make sure initvals are within the bounds
lowind = (initvals<lb);
upind = (initvals>ub);

for i = 1:length(initvals)
    
    if lowind(i) ~= 0
        fprintf('paramter %i is less than it''s lower bound\n',i)
    end
    
end

for i = 1:length(initvals)

    if upind(i) ~= 0
        fprintf('paramter %i is larger than it''s upper bound\n',i)
    end
    
end

if norm(abs(lowind)) ~= 0 || norm(abs(upind)) ~= 0
    error('there are initial values that are outside of the bounds')
end


% extract exectuable name from fortran_code
if brutus
    execname = 'hc_v02_J';
else
    fexe = dir('fortran_code/*.exe');
        
    if isempty(fexe) 
        error('custom error:: no *.exe file found')
    end
    
    start_ind = regexp(fexe.name,'.exe','once'); % cut the extension
    execname = fexe.name(1:start_ind-1); % use short names for brutus
end

% initialize log-file
% STOP-FLAG ------------------ DO NOT DELETE ------------------------------

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
    options.PlotFcns = @psplotbestf;
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
outrowvec = [optimoutput.params'.*initvals'...
    optimoutput.resnorm optimoutput.flag];
fprintf(fid,'%s\n',headerfunc(1:numpars));
fprintf(fid,repmat('%f\t',1,numcols),outrowvec);
fclose(fid);