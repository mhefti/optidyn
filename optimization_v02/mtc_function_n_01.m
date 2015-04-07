%% mass transfer scaling
% contains a number of (empirical) models specified in'type' for the mass
% transfer scaling function in terms of the relative humidity
% outputs a <figure_name_type>.png figure where 'type' is the model

function mtc_function_n_01

% INPUTS

% get from optim_sim
fpath = strcat('output_optim',filesep,'worker_01',filesep);
pfun = @(fname) strcat(fpath,fname);

% get isotherm type and parameters
fid = fopen(pfun('isotherm.dat'),'r+');
isotype = textscan(fid,'%s','delimiter','\t','HeaderLines',0);
isotype = isotype{1}(1);
isotype = sprintf('%s',isotype{:});
frewind(fid)
fcontent = textscan(fid,'%f','delimiter','\t','HeaderLines',1);
fclose(fid);

if ~isempty(regexp('''Sips_Sips''',isotype,'once'))
    parw = fcontent{1}(1:6);
elseif ~isempty(regexp('''Do_modified''',isotype,'once'))
    parw = fcontent{1}(1:5);
else
    error('unknown isotherm type')
end

% get mtc model type
fid = fopen(pfun('settings.dat'),'r+');
% ensure every line has an 'end of line character'; enable eol in Notepad++
% to check
clear fcontent
fcontent = textscan(fid,[repmat('%s',1,10) '%*[^\n]']);
fclose(fid);
mtcmodelid = 16; % retrieve from expinfo in future versions
typee = char(fcontent{1}(mtcmodelid));

% get mtc0
formatSpec = '%f%[^\n\r]';
fid = fopen(pfun('conditions.dat'),'r+');
clear fcontent
fcontent = textscan(fid, formatSpec, 'Delimiter', '\t',  'ReturnOnError', false);
fclose(fid);
mtc0id = 16;
mtc0 = fcontent{1}(mtc0id);

% get mtc-model parameters
formatSpec = '%s%[^\n\r]';
fid = fopen(pfun('fitting.dat'),'r+');
clear fcontent
fcontent = textscan(fid, formatSpec, 'Delimiter','','ReturnOnError',false);
fclose(fid);
nmtc = str2double(cell2mat(fcontent{1}(1)));
mtcscale = zeros(nmtc,1);
j = 0;

for i=1:nmtc
    mtcscale(i) = str2double(cell2mat(fcontent{1}(3 + j)));
    j = j + 2;
end

% END OF INPUTS

fprintf('\n')
fprintf('------------------------------------------\n')
fprintf('mtc-model:                     %s\n',typee)
fprintf('# of mtc parameters            %i\n',nmtc)
fprintf('mtc0:                          %f\n',mtc0)
for i=1:nmtc
fprintf('mtc_scale_%i                    %f\n',i,mtcscale(i))
end
fprintf('isotherm-type                 %s\n',isotype)
fprintf('------------------------------------------\n')

nstar = linspace(0,1,1e3);
nofrv = @(rv) rv;

switch typee
    
    case '''monotonic'''
        
        scaler = @(rr)  1 - mtcscale(1)*nofrv(rr).^mtcscale(2)./...
            (1 + mtcscale(1)*nofrv(rr).^mtcscale(2));
        
    case '''monoplus'''
        
        scaler = @(rr)  1 + mtcscale(1)*mtcscale(2)*nofrv(rr).^mtcscale(3)./...
            (1 + mtcscale(2)*nofrv(rr).^mtcscale(3));
        
    case '''exponential'''
        
        scaler = @(rr) mtcscale(1)*exp(-mtcscale(2)*nofrv(rr)) + 1;
        
    case '''custom'''
        
        scaler = @(rr) mtcscale(1)*exp(-mtcscale(2)*nofrv(rr))+ mtcscale(3)*mtcscale(4)*nofrv(rr).^mtcscale(5)./...
            (1 + mtcscale(4)*nofrv(rr).^mtcscale(5)) + 1;
end

mtcdat.x = nofrv(nstar)*isotherm(char(isotype),parw,1);
mtcdat.mtc = mtc0*scaler(nstar);

figure(1)
clf
hold on
plot(mtcdat.x,mtcdat.mtc,'k')
xlabel('n [mol/kg]')
ylabel('mtc [s^{-1}]')
export_fig(strcat('mtc_func_',typee,'.png'))
set(gca,'yscale','log')
export_fig(strcat('mtc_func_',typee,'_logy.png'))

mtcdat.mtc0 = mtc0; 
mtcdat.mtcscale = mtcscale; 
mtcdat.isotype = isotype;
mtcdat.isopar = parw; 
mtcdat.mtcmodel = typee; 
mtcdat.rh = linspace(0,1,5e2); 
mtcdat.n = isotherm(isotype,parw,mtcdat.rh);

figure(2) 
clf 
hold on
plot(mtcdat.rh,mtcdat.n); 
xlabel('rh [-]')
ylabel('n [mol/kg]')
export_fig 'isotherm.png'

save('mtcdat','mtcdat')


function [xi,term1,term2] = isotherm(isotype,parw,x)
    
    switch isotype
        
        case '''Sips_Sips'''
            S0 = parw(1); k = parw(2); m0 = parw(3);
            Smu = parw(4); Kmu = parw(5); m = parw(6);
            
        case '''Do_modified'''
            S0 = parw(1); k = parw(2); m0 = 1;
            Smu = parw(3); Kmu = parw(4); m = parw(5);
            
    end
    
    term1 = S0*k*x.^m0./(1 + k*x.^m0);
    term2 = Smu*Kmu*x.^m./(1 + Kmu*x.^m);
    xi = term1 + term2;
    
end

end









