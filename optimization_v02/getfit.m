function output = getfit(hand,initvals,x,y,lb,ub,opt,execname,alg,parall,brut,expinfo,slaves)

% check inputs
if nargin >= 10
    
    if strcmp(parall,'parallel')
        parall = true;
    else
        parall = false;
    end
    
    if strcmp(brut,'brutus')
        brutus = true;
    else
        brutus = false;
    end
    
else
    parall = false;
    brutus = false;
end

dirfcn = @(ii) sprintf('worker_%2.2d',ii);
% -------------------------------------------------------------------------
% prepare for parallel jobs    
if parall
        
    if brutus 
        
        if ~expinfo.eulermode % runs on Brutus
            poolobj = parcluster('BrutusLSF8h');
        else % runs on Euler
            poolobj = parcluster('EulerLSF8h');
        end
    
        sarg = sprintf('-W %s:00 -R "rusage[mem=1000]"',num2str(expinfo.time));
        poolobj.SubmitArguments = sarg;
        jobid = getenv('LSB_JOBID');
        mkdir(jobid);
        poolobj.JobStorageLocation = jobid;
        poolobj.NumWorkers = slaves;
    
        % sometimes the pool fails to open
        while (matlabpool('size') <= 0)
            
            try
                matlabpool(poolobj,poolobj.NumWorkers)
            catch
            end
            
        end
        
    else % parallel on local machine
        
        if matlabpool('size') ~= 0
            matlabpool('close');
        end
        
        poolobj = parcluster;
        poolobj.NumWorkers = 2;
        matlabpool(poolobj)
        
        for i=1:poolobj.NumWorkers
            
            if exist(dirfcn(i),'dir')
                
                try
                    rehash();
                    taskkill = sprintf('taskkill /F /IM %s.exe',execname);
                    system(taskkill)
                    rmdir(dirfcn(i),'s');
                catch
                end
                
            end
            
        end
        
    end % / brutus
    
    for i=1:poolobj.NumWorkers
        
        if exist(dirfcn(i),'dir') % leave for brutus
            rehash();
            % rmdir(dirfcn(i),'s');
        end
        
        copyfile('prms/*.dat',dirfcn(i));
        copyfile('dlmcell.m',dirfcn(i));
        
        if expinfo.fitiso % copy contents of isothermfit if enabled
            copyfile('isotherm_fit/*',dirfcn(i));
        end
        
        % copy executable to worker_iD folder
        if brutus
            copyfile(sprintf('fortran_code/%s',execname),dirfcn(i));
            % rename executable, attach iD at end of filename
            movefile(strcat(dirfcn(i),filesep,execname),...
                strcat(dirfcn(i),filesep,sprintf('%s_%2.2d',execname,i)));
        else
            copyfile(sprintf('fortran_code/%s.exe',execname),dirfcn(i));
            % rename executable, attach iD at end of filename
            movefile(strcat(dirfcn(i),filesep,execname,'.exe'),...
                strcat(dirfcn(i),filesep,sprintf('%s_%2.2d.exe',execname,i)));
        end
        
   end
       
else    % if not parallel, the action happens in directory worker_01
        
    if exist(dirfcn(1),'dir')

        try
            rehash();
            % rmdir(dirfcn(i),'s');
            
            if brutus 
                taskkill = sprintf('pkill %s_%2.2d.exe',execname,1); 
            else
                taskkill = sprintf('taskkill /F /IM %s_%2.2d.exe',execname,1);
            end
            
            system(taskkill);
        catch
        end
        
    end
    
    copyfile('prms/*.dat',dirfcn(1));
    copyfile('dlmcell.m',dirfcn(1));
    
    if expinfo.fitiso % copy contents of isothermfit if enabled
        copyfile('isotherm_fit/*',dirfcn(1));
    end
    
    if brutus % no .exe extension
        copyfile(sprintf('fortran_code/%s',execname),dirfcn(1));
        % rename executable, attach iD at end of filename
        movefile(strcat(dirfcn(1),filesep,execname),...
            strcat(dirfcn(1),filesep,sprintf('%s_%2.2d',execname,1)));
    else
        copyfile(sprintf('fortran_code/%s.exe',execname),dirfcn(1));
        % rename executable, attach iD at end of filename
        movefile(strcat(dirfcn(1),filesep,execname,'.exe'),...
            strcat(dirfcn(1),filesep,sprintf('%s_%2.2d.exe',execname,1)));
    end
    
    
end % / parall


% define the objective function

optimfcn = @(param) optimfcn_sa(param,y,hand,execname,expinfo,brutus,initvals);

% for normalization, input ones to solver
initvals_optim = ones(size(initvals));
lb = lb./initvals;
ub = ub./initvals;

switch alg 
    
    case 'lsq'
        opt = optimset('MaxFunEvals',1e4);
        [output.params,output.resnorm,output.residual] = ...
            lsqcurvefit(hand,initvals_optim,x,y,lb,ub,opt);

    case 'MultiStart'
        problem = createOptimProblem('fmincon','objective',optimfcn,...
            'x0',initvals_optim,'lb',lb,'ub',ub,'options',opt);
        
        if parall
            ms = MultiStart('UseParallel','always');
            [output.params,output.resnorm,output.flag,...
                output.out,allmins] = run(ms,problem,expinfo.numMSpoints);
            matlabpool close
        else
            ms = MultiStart('UseParallel','never');
            [output.params,output.resnorm,output.flag,...
                output.out,allmins] = run(ms,problem,100);
        end
        
    case 'GlobalSearch'
        problem = createOptimProblem('fmincon','objective',optimfcn,...
            'x0',initvals_optim,'lb',lb,'ub',ub,'options',opt);
        gs = GlobalSearch;
        [output.params,output.resnorm,output.flag,...
            output.out,allmins] = run(gs,problem);
                
    case 'fminsearchbnd'
        [output.params,output.resnorm,output.flag,output.out] = ...
            fminsearchbnd(optimfcn,initvals_optim,lb,ub,opt);        
        
    case 'patternsearch'
        
        if parall
            opt.UseParallel = 'always';
            opt.Vectorized = 'off';
            opt.CompletePoll = 'on';
        else
            opt.UseParallel = 'never';
        end
        
        [output.params,output.resnorm,output.flag,output.out] = ...
            patternsearch(optimfcn,initvals_optim,[],[],[],[],lb,ub,[],opt);
        
end