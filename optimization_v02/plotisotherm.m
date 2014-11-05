function plotisotherm(params,leginfo)
% tiny function to return a figure with the same # of isotherms as # of 
% parent fields contained in the input structure params

if nargin <1
    params = [];
end

if nargin < 2
    leginfo =[];
end

datasource = 'D:\Users\mhefti\Documents\Isotherm Models\Activated Carbons\AP3-60_new_sensor_balance\';
acdata = load(strcat(datasource,'acdata'));
params_exp_fit = load(strcat(datasource,'param_DoLang_ads'));
exp_fit = params_exp_fit.param_DoLang_ads;
Tshift = 0;
acdata.ads.RH = Pvap_H2O(acdata.ads.Tdew,'Pa')./Pvap_H2O(acdata.ads.T + Tshift,'Pa');
% exp_fit = [1114	0.0031619	21.141	33.602	9.5747];


numseries = length(params);
RHvec = linspace(0,1,1e3);


figure(1)
clf
box on; grid on; hold on
set(gca,'LineWidth',1.1,'TickLength',[0.008 0.008],'FontSize',14)
xlabel('relative humidity [%]')
ylabel('adsorbed amount [mol/kg]')
xlim([0 1])
ylim([0 25])

% adsorption data
nadspoints = numel(acdata.ads.n);
rr_fit = nadspoints - 13:nadspoints;
rr_rest = 1:nadspoints-13;
fcol = [128,0,0]/255;
rcol = [128,128,128]/255;
plot(acdata.ads.RH(rr_fit),acdata.ads.n(rr_fit),'o','Color',fcol,'MarkerFaceColor',fcol)
plot(acdata.ads.RH(rr_rest),acdata.ads.n(rr_rest),'o','Color',rcol,'MarkerFaceColor',rcol)
% desorption data
% plot(acdata.des.RH,acdata.des.n,'o k','MarkerFaceColor','w')

% fits to experimental data
plot(RHvec,model_DoLang(exp_fit,RHvec),'k','LineWidth',2)


colmat = [254,217,118; 254,178,76;253,141,60; 240,59,32; 189,0,38];
colmat = colmat(:,end:-1:1);
colmat = colmat/255;

for i=1:numseries
	plot(RHvec,model_DoLang(params(i).series,RHvec),'Color',colmat(i,:),'LineWidth',2)
end



if isempty(leginfo) ~= true
    [rw,cl] = size(leginfo);
    
    if cl > rw 
        leginfo = leginfo';
    end
    
end

lcell = [{'exp #: 140828';'static rest';'fit to 140828'}; leginfo];
legh = legend(lcell,'Location','North');
set(legh,'FontSize',10)
set(legh,'Box','off')
print('-depsc','isotherms.eps')


function xi = model_DoLang(parw,x)

S0 = parw(1); k = parw(2);
Smu = parw(3); Kmu = parw(4); m = parw(5);
xi = S0*k*x./(1 + k*x)...
    + Smu*Kmu*x.^m./(1 + Kmu*x.^m);
end


end