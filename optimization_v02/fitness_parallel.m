function outdata = fitness_parallel(params,execname,expinfo,brutus,initvals)
% computes output of the model
% timeout function built in, set timeout below

task = getCurrentTask();

if ~isempty(task)
    task_curr = task.ID;
else
    task_curr = 1;
end

% name functions
sfield = @(q,ii) sprintf('%s_%s_%s',q,num2str(expinfo.list(ii)),expinfo.modes(ii));

% if ~expinfo.fiteq % detailed model is fitted
fprintf('----------------------------------------------------\n')
fprintf('current worker ID:  %d\n',task_curr)
wd_curr = sprintf('worker_%2.2d',task_curr);
fprintf('current parameters evaluated: %f\n',params)
% change directory
cd(wd_curr) 
% remove the dump file
succ_dump = 'dump_succ.log';
%end


if exist(succ_dump,'file')
    delete(succ_dump);
end

%% parameter handling

% read/rewrite the files % replace just w vector and add to expinfo!

fid = fopen('isotherm.dat','r+');
fcontent = textscan(fid,'%f','delimiter','\t','HeaderLines',1);

% parse the isotherm parameters to be fitted
if expinfo.isoflex
   
    for i=1:numel(expinfo.fitisolocs)
        fcontent{1}(expinfo.fitisolocs(i)) = params(i);    
    end

else
    fcontent{1}(1:expinfo.num_isopar) = params(1:expinfo.num_isopar);
end

fcontent = [fcontent{1}(:)];
frewind(fid)
fprintf(fid,'%s\n',strcat('''',expinfo.iso_type,'''')); % jap, 4 apostrophes needed
fprintf(fid,'%f\n',fcontent);
fprintf(fid,'%s\n','''Inert''');
fclose(fid);

% parameter fitting
if expinfo.fitiso 
    outdata.iso_pars = isotherm_fit(params(1));
end

% heat of adsorption
ind_loc = expinfo.num_pars;
Hads = [];
if expinfo.fit_isothermal == true && expinfo.fit_Hads == true
    Hads = 0;
    warning('Hads was set to zero because isothermal fitting enabled')
elseif expinfo.fit_Hads
    Hads = params(ind_loc); % always the last entry if enabled
    ind_loc = expinfo.num_pars - 1; % move one up in the fitting vector
end

if ~isempty(Hads)
    fid = fopen('parameter1.dat','r+');
    fcontent = textscan(fid,'%s%s %*[^\n]');
    fclose(fid);
    fcontent{1}(expinfo.Hads_id) = num2cell(-Hads*1e3);
    dlmcell('parameter1.dat',[fcontent{1} fcontent{2}]);
end

% heat transfer coefficient, mass transfer coefficient
htc = [];
if expinfo.fit_isothermal == true && expinfo.fit_htc == true
    htc = 1e6;
    warning('htc was set to 1e6 because isothermal fitting enabled')
elseif expinfo.fit_htc == true 
    htc = params(ind_loc);
    ind_loc = ind_loc - 1; % move one up in the fitting vector
end

mtc = [];
if expinfo.fit_mtc
    mtc = params(ind_loc);
    ind_loc = ind_loc - 1;
end
    
% update conditions for fitting
newconditions = expinfo.conditions;

if ~isempty(htc)
    newconditions(expinfo.htc_id,:) = htc; % htc
end

if ~isempty(mtc)
    newconditions(expinfo.mtc_id,:) = mtc; % htc
end

if isempty(htc) == false || isempty(mtc) == false % conditions rewritten
    dlmwrite('conditions.dat',newconditions,'\t')
end

%% computation of the model

% timeout function
timeout = expinfo.time_out; % s
tstart = tic;
% filenames
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
