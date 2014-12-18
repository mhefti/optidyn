% read fitting log
optimpars = readfittinglog;

% create the folder structure
mkdir('output_optim')
mkdir('output_optim/worker_01')

% get the .exe and prms files
copyfile('prms/*.dat','output_optim/worker_01/')
copyfile('fitness_parallel.m','output_optim/')
copyfile('dlmcell.m','output_optim/worker_01/')

% prepare for plotting
gppath = 'D:\Users\mhefti\Documents\Projects\dynamic_experiments_modelling\20141016_hysteresis_clean_v00_jacobian_v02\prms\plot_simulation_experiment_mult_comparison.plt';
copyfile(gppath,'output_optim/worker_01')

% extract part of the optimization_dynamic_v02
copyfile('optimization_dynamic_v02.m','temp.txt')
content = fileread('temp.txt');
[~,ind_t] = regexp(content,'STOP-FLAG');

% replace brutmode = 'brutus' to brutmode = 'nobrutus';
expr = sprintf('%s %s %s','brutmode','=','''brutus''');
[ind_s, ind_e] = regexp(content,expr);

if ~isempty(ind_s)
    [ind_s] = regexp(content,'''brutus''');
    ind_s = ind_s(1);
    numchar = length('brutus') + 2;
    content(ind_s:ind_s+numchar+1) = '''nobrutus''';
end

dlmwrite('temp.m',content(1:ind_t),'')
run('temp.m')
% save some workspace; expinfo available
save('output_optim/temp_space')
copyfile(strcat('fortran_code/',execname,'.exe'),...
    sprintf('output_optim/worker_01/%s_%2.2d%s',execname,1,'.exe'));

cd 'output_optim/'
load('temp_space')
[outdata] = fitness_parallel(optimpars(1:end-1),execname,expinfo,false,[]);
fprintf('the resnorm of the optimum was                 %f\n',optimpars(end))  

% plot the results
cd 'worker_01/'
nn = dir('*.plt');
comm = sprintf('gnuplot %s',nn.name);
system(comm)
cd ../..

delete('temp.m')