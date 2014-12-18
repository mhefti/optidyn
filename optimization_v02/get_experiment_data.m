function outstruc = get_experiment_data(dateval,brutmode,timecut,mode_in,mtc_0) 
% reads the experiment output at the date = dateval from the T drive
% input: dateval is the date of form YYMMDD as number 
% outputs: - outstruc with time, concentration and temperature (at z = 0.5) 
%          - writes conditions file for the fortran simulation
if nargin == 2
    timecut = [];
end
%% set up paths & names
if strcmp(brutmode,'brutus')
    
    destpath = strcat('..',filesep,'experiments_by_date',filesep);
    
    if exist(destpath,'dir')
        destpath = destpath; % do nothing
    elseif exist(strcat('..',filesep,destpath),'dir')
        destpath = strcat('..',filesep,destpath);
    else
        error('cannot find experiments_by_date directory')
    end
        
else
    destpath = '\\d.ethz.ch\dfs\groups\mavt\spl\Temp\adsorption\SmallColumn\eval\experiments\experiments_by_date\';
end

simpath = '';
datein = dateval;
datestr = num2str(datein);
destpath = strcat(destpath,datestr,filesep);

if strcmp(mode_in,'a') 
    mode = 'adsorption';
    adsmode = true;
else
    mode = 'desorption';
    adsmode = false;
end
file_LV = @(type) strcat(destpath,'outputLV_',type,'_',datestr,'_',mode,'.dat');
file_dP = @(type) strcat(destpath,type,'_',datestr,'_',mode,'.dat');
%% outstruc
outstruc = read_file(file_LV('column'),1,'column');
% find the real starting time, assume is at minimum of yH2O:
[~,Imin] = min(outstruc.yH2O);
[~,Imax] = max(outstruc.yH2O);

% if timecut not specified take the experiments from t_start to end
if isempty(timecut)
    timecut = max(outstruc.time);
end

dt = outstruc.time(end) - outstruc.time(end-1);
% find the location of timecut in outstruc.time
[tmax,~] = find(abs((outstruc.time - timecut)) < dt/exp(1));

set_zero = 1:Imin-1;
set_max = 1:Imax-1;
therange = 1:tmax;

%% write files
% relevant file names
data_LV_Tcol = read_file(file_LV('Tcol'),1,'Tcol');
data_LV_means = read_file(file_LV('means'),1,'means');
data_dP = read_file(file_dP('dP'),1,'dP');
% conditions --------------------------------------------------------------
i = 1;
% time 
time = outstruc.time(therange); condmat(i,1) = max(time); i = i + 1;
% therange = 1:max(time);
tshift = 0; condmat(i,1) = tshift; i = i + 1;
% flow
flow = data_dP.flow(1)*1e-6; condmat(i,1) = flow; i = i + 1;
% pressure 
Pfeed = data_dP.pmeasure(1); condmat(i,1) = Pfeed; i = i + 1;
P0 = Pfeed; condmat(i,1) = P0; i = i + 1;
Pamb = Pfeed; condmat(i,1) = Pamb; i = i + 1;
% temperature
Tfeed = mean(data_LV_Tcol.T(end-10:end)) + 273.15; condmat(i,1) = Tfeed; i = i + 1;
T0 = Tfeed; condmat(i,1) = T0; i = i + 1;
Tamb0 = Tfeed; condmat(i,1) = Tamb0; i = i + 1;
Tambend = Tfeed; condmat(i,1) = Tambend; i = i + 1;
% theat
theat = 100.0; condmat(i,1) = theat; i = i + 1;
% yfeed
yfH2O = data_LV_means.yfeed; condmat(i,1) = yfH2O; i = i + 1;
yfHe = 1.0 - yfH2O; condmat(i,1) = yfHe; i = i + 1;
% y0
if adsmode 
    y0H2O = min(outstruc.yH2O); 
else
    y0H2O = max(outstruc.yH2O);
end

condmat(i,1) = y0H2O; i = i + 1;
y0He = 1 - y0H2O; condmat(i,1) = y0He; i = i + 1;
% mtc 
kH2O = mtc_0; condmat(i,1) = kH2O; i = i + 1;
kHe = 1; condmat(i,1) = kHe; i = i + 1;
% hL
hL = 220; condmat(i,1) = hL; i = i + 1;
% exp #
exnum = 1; condmat(i,1) = exnum;
% parameter2 --------------------------------------------------------------
i = 1; 
UW = 200; p2mat(i,1) = UW; i = i + 1;
UWA = 10; p2mat(i,1) = UWA; i = i + 1;
H2Oid = '''H2O'''; %p2mat(i,1) = H2Oid; i = i + 1;
Heid = '''He'''; %p2mat(i,1) = Heid; i = i + 1;
dlmwrite(strcat(simpath,'prms/parameter2.dat'),p2mat)
dlmwrite(strcat(simpath,'prms/parameter2.dat'),H2Oid,'delimiter','','-append')
dlmwrite(strcat(simpath,'prms/parameter2.dat'),Heid,'delimiter','','-append')

% add the conditions to the outstruc
if adsmode 
    outstruc.yH2O(set_zero) = y0H2O*ones(max(set_zero),1);
else
    outstruc.yH2O(set_max) = outstruc.yH2O(Imax)*ones(set_max,1);
end

outstruc.t = outstruc.time(therange);
outstruc.yH2O = outstruc.yH2O(therange);
outstruc.T05 = data_LV_Tcol.T(therange);

if  strcmp(brutmode,'brutus') == false
    figure('Name',sprintf('%s %s',datestr,mode))
    clf; grid on;
    xlabel('time [min]')
    ylabel('yH2O [-]')
    plot(therange,outstruc.yH2O,'k','LineWidth',1.5)
end


outstruc.condmat = condmat;