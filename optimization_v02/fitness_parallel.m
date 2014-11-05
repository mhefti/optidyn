function outdata = fitness_parallel(params,execname,expinfo,brutus,initvals)
% computes output of the model
% timeout function built in, set timeout below

task = getCurrentTask();

if ~isempty(task)
    task_curr = task.ID;
else
    task_curr = 1;
end

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

% read/rewrite the files % replace just w vector and add to expinfo!
fid = fopen('isotherm.dat','r+');
fcontent = textscan(fid,'%f','delimiter','\t','HeaderLines',1);
fcontent{1}(1:5) = params(1:5);
fcontent = [fcontent{1}(:)];
frewind(fid)
fprintf(fid,'%s\n','''Do_modified''');
fprintf(fid,'%f\n',fcontent);
fprintf(fid,'%s\n','''Inert''');
fclose(fid);

% parameter fitting
if expinfo.fitiso 
    outdata.iso_pars = isotherm_fit(params(1));
end

% update conditions for fitting
newconditions = expinfo.conditions;
exptimes = newconditions(1,:);
newconditions(expinfo.fparamids(1),:) = params(6); % mtc
% newconditions(expinfo.fparamids(2),:) = params(7); % mtc
dlmwrite('conditions.dat',newconditions,'\t')

% for heat of adsorption
fid = fopen('parameter1.dat','r+');
fcontent = textscan(fid,'%s%s %*[^\n]');
fclose(fid);
fcontent{1}(expinfo.fparamids(2)) = num2cell(-params(7));
dlmcell('parameter1.dat',[fcontent{1} fcontent{2}]);

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
        fprintf('-------------------TIMEOUT----------------------------\n')
        timeup = true;
        running = false;
    end
    
end

fclose(fid);

%% generate the function output
sfield = @(q,ii) sprintf('%s_%s_%s',q,num2str(expinfo.list(ii)),expinfo.modes(ii));

if timeup || noerr_flag ~= true
    
    for i = 1:expinfo.numexp
        
        if expinfo.fity
            outdata.(sfield('yH2O',i)) = pi*ones(expinfo.numpoints(i),1); 
        end
        
        if expinfo.fitT
            outdata.(sfield('T05',i)) = pi*ones(expinfo.numpoints(i),1); 
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
