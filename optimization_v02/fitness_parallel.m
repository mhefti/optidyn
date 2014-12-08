function outdata = fitness_parallel(params,execname,expinfo,brutus,initvals)
% computes output of the model
% timeout function built in, set timeout below

task = getCurrentTask();

if ~isempty(task)
    task_curr = task.ID;
else
    task_curr = 1;
end

% if ~expinfo.fiteq % detailed model is fitted
fprintf('----------------------------------------------------\n')
fprintf('current worker ID:  %d\n',task_curr)
wd_curr = sprintf('worker_%2.2d',task_curr);
fprintf('current parameters evaluated: %f\n',params)

% change directory
cd(wd_curr) 

% remove the dump file
succ_dump = 'dump_succ.log';

if exist(succ_dump,'file')
    delete(succ_dump);
end

%% parameter handling
% *** NOTE: keep order of the variables, as specified in the file 
% optimization_dynamic_v0X
% -------------------------------------------------------------------------
% [     isopar(1)                   % 1
%          .                        % .
%          .                        % .
%       isopar(numisopar)           % numisopar
%       mtc                         % numisopar + 1
%       mtcmodelpar_1               % numisopar + 2
%       .                           % .
%       mtcmodelpar_2               % numisopar + nummtcpar
%       htc                         % numisopar + nummtcpar + 1
%       Hads                ]       % numpar
% -------------------------------------------------------------------------

% isotherm parameters
% -------------------------------------------------------------------------
fid = fopen('isotherm.dat','r+');
fcontent = textscan(fid,'%f','delimiter','\t','HeaderLines',1);

% parse the isotherm parameters to be fitted
if expinfo.fitiso % dynamic isotherm parameter fitting to static data
    
    if expinfo.isoflex
        % get the parameters that are fixed from isotherm.dat, i.e. all others
        % than specified in expinfo.fitisolocs
        fix_idx = setdiff(1:expinfo.num_isopar,expinfo.fitisolocs);
        fixisopars = fcontent{1}(fix_idx);
        outdata.iso_pars = isotherm_fit(params(1),expinfo,[fixisopars],fix_idx);
        
        for i=1:numel(expinfo.fitisolocs)
            fcontent{1}(expinfo.fitisolocs(i)) = outdata.iso_pars(i);
        end
        
    else
        outdata.iso_pars = isotherm_fit(params(1),expinfo,[],[]);
        fcontent{1}(1:numel(outdata.iso_pars)) = outdata.iso_pars;
    end
    
else % isotherm fitting to dynamic experiments directly

    if expinfo.isoflex
   
        for i=1:numel(expinfo.fitisolocs)
            fcontent{1}(expinfo.fitisolocs(i)) = params(i);    
        end

    else
        fcontent{1}(1:expinfo.num_isopar) = params(1:expinfo.num_isopar);
    end 
    
end

% write the contents of isotherm.dat
fcontent = [fcontent{1}(:)];
frewind(fid)
fprintf(fid,'%s\n',strcat('''',expinfo.iso_type,'''')); % jap, 4 apostrophes needed
fprintf(fid,'%f\n',fcontent);
fprintf(fid,'%s\n','''Inert''');
fclose(fid);

% all other parameters to be fitted; 
ind_loc = expinfo.num_pars; % start at the end of the par vector

% **NOTE: keep the order of the parameters

% heat of adsorption 
% -------------------------------------------------------------------------
Hads = [];
if expinfo.fit_isothermal == true && expinfo.fit_Hads == true
    Hads = 0;
    warning('Hads was set to zero because isothermal fitting enabled')
elseif expinfo.fit_Hads
    Hads = params(ind_loc); % always the last entry if enabled
    ind_loc = ind_loc - 1; % move one up in the fitting vector
end

if ~isempty(Hads)
    fid = fopen('parameter1.dat','r+');
    fcontent = textscan(fid,'%s%s %*[^\n]');
    fclose(fid);
    fcontent{1}(expinfo.Hads_id) = num2cell(-Hads*1e3);
    dlmcell('parameter1.dat',[fcontent{1} fcontent{2}]);
end

% heat transfer coefficient, mass transfer coefficient
% -------------------------------------------------------------------------
htc = [];
if expinfo.fit_isothermal == true && expinfo.fit_htc == true
    htc = 1e6;
    warning('htc was set to 1e6 because isothermal fitting enabled')
elseif expinfo.fit_htc == true 
    htc = params(ind_loc);
    ind_loc = ind_loc - 1; % move one up in the fitting vector
end

% mass transfer
% -------------------------------------------------------------------------
mtc = [];
if expinfo.fit_mtc
    mtc = params(ind_loc);
    ind_loc = ind_loc - 1;
end

fid = fopen('settings.dat','r+');
fcontent = textscan(fid,[repmat('%s',1,8) '%*[^\n]']);
fclose(fid);
fcontent{1}(expinfo.mtcmodelid) = {strcat('''',expinfo.mtcmodel,'''')};
dlmcell('settings.dat',[fcontent{1} fcontent{2:end}]);

% mass transfer model - write in fitting.dat
if expinfo.fitmtcmodel
    % write the mtcmodel parameters in fitting.dat 
    mtcstr = @(ii) sprintf('''mtcscale_%i''',ii);
    mtccell = cell(2*expinfo.nummtcpar,1);
    j = 2*expinfo.nummtcpar;
    
    % needs to be written 'backwards'
    for i = expinfo.nummtcpar:-1:1  
        mtccell{j} = params(ind_loc);
        mtccell{j - 1} = mtcstr(i);
        j = j - 2;
        ind_loc = ind_loc - 1;
    end

    % write the parameters in fitting.dat
    mtccell = [num2str(expinfo.nummtcpar); mtccell];
    dlmcell('fitting.dat',mtccell);  
end
      
% update conditions for fitting
newconditions = expinfo.conditions;

if ~isempty(htc)
    newconditions(expinfo.htc_id,:) = htc;
end

if ~isempty(mtc)
    newconditions(expinfo.mtc_id,:) = mtc;
end

if isempty(htc) == false || isempty(mtc) == false % conditions rewritten
    dlmwrite('conditions.dat',newconditions,'\t')
end

%% computation of the model
% timeout function
timeout = expinfo.time_out; % s
tstart = tic;
logfname = 'simlog.log';

if brutus
    % no .exe extension in brutus environment
    execname_curr = sprintf('%s_%2.2d',execname,task_curr);
    system(sprintf('./%s>%s &',execname_curr,logfname));
else
    execname_curr = sprintf('%s_%2.2d.exe',execname,task_curr);
    system(sprintf('start /B %s>%s & cmd.exe /C',execname_curr,logfname));%2>&1    
end

running = true;
timeup = false;

pause(0.1);
fid = fopen('simlog.log','r');

while running == true && exist(succ_dump,'file') == false
    runtime = toc(tstart);
    simlog = textscan(fid,'%s');
    frewind(fid);
    % check if ERROR was detected in output
    noerr_flag = isempty(cell2mat(regexp(simlog{:},'ERROR','once')));
    
    if noerr_flag ~= true
        running = false;
        fprintf('--------------------ERROR-----------------------------\n')
    end
    
    if runtime > timeout 
        
        if brutus
            taskkill = sprintf('pkill %s',execname_curr);
        else
            taskkill = sprintf('taskkill /F /IM %s',execname_curr); 
        end
        
        system(taskkill);
        fprintf('-------------------TIMEOUT--------------------------\n')
        timeup = true;
        running = false;
    end
    
end

fclose(fid);

%% generate the function output
sfield = @(q,ii) sprintf('%s_%s_%s',q,num2str(expinfo.list(ii)),expinfo.modes(ii));

if timeup || noerr_flag ~= true
    
    for i = 1:expinfo.numexp
        outdata.(sfield('time',i)) = linspace(1,expinfo.numpoints(i),expinfo.numpoints(i));
        
        if expinfo.fity
            outdata.(sfield('yH2O',i)) = linspace(pi,pi*pi,expinfo.numpoints(i))'; 
        end
        
        if expinfo.fitT
            outdata.(sfield('T05',i)) = linspace(pi,pi*pi,expinfo.numpoints(i))'; 
        end
        
    end
    
else
    % retrieve data from simulation output
    expfolder = @(ii) sprintf('experiment_%2.2d',ii);
    
    for i = 1:expinfo.numexp
        
        if expinfo.fity
            % reading yH2O
            fid = fopen(strcat(expfolder(i),'/exitprofile.txt'),'r');
            simdata = textscan(fid,repmat('%f',1,9),'delimiter','\t',...
                'HeaderLines',1);
            fclose(fid);
            outdata.(sfield('yH2O',i)) = simdata{2}(:);
            outdata.(sfield('time',i)) = simdata{1}(:);
        end
        
        if expinfo.fitT
            % reading T05
            fid = fopen(strcat(expfolder(i),'/temperatures.txt'),'r');
            simdata = textscan(fid,repmat('%f',1,2),'delimiter','\t',...
                'HeaderLines',1);
            fclose(fid);
            outdata.(sfield('T05',i)) = simdata{2}(:);
            outdata.(sfield('time',i)) = simdata{1}(:);
        end
        
    end
    
end

cd ..
