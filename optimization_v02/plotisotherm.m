function plotisotherm(params,leginfo,type)
% tiny function to return a figure with the same # of isotherms as # of 
% parent fields contained in the input structure params

if nargin <1
    params = [];
end

if nargin < 2
    leginfo =[];
end

if nargin < 3
    type = 'sips_sips';
end

datasource = 'D:\Users\mhefti\Documents\Isotherm Models\Activated Carbons\AP3-60_new_sensor_balance\';
acdata = load(strcat(datasource,'acdata'));
params_exp_fit = load(strcat(datasource,'param_DoLang_ads'));
exp_fit = params_exp_fit.param_DoLang_ads;
Tshift = 0;
acdata.ads.RH = Pvap_H2O(acdata.ads.Tdew,'Pa')./Pvap_H2O(acdata.ads.T + Tshift,'Pa');
acdata.des.RH = Pvap_H2O(acdata.des.Tdew,'Pa')./Pvap_H2O(acdata.des.T + Tshift,'Pa');
% exp_fit = [1114	0.0031619	21.141	33.602	9.5747];


numseries = length(params);
RHvec = linspace(0,1,1e3);


figure(1)
clf
box on; grid on; hold on
set(gca,'LineWidth',1.1,'TickLength',[0.008 0.008],'FontSize',14)
xlabel('relative humidity [-]')
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
plot(acdata.des.RH,acdata.des.n,'o k','MarkerFaceColor','w')

% fits to experimental data
plot(RHvec,model_DoLang(exp_fit,RHvec),'k','LineWidth',2)



% source: D:\Users\mhefti\Documents\Projects\dynamic_exp_mod\dynamic_exp_fit\
% optimization_v02\EULER_by_date\141219_dispatch_jobs\141215_exp140714_
%140718_SS_all_mtc_custom_fity_MLE_120\output_optim\worker_01

opt_fit = [3084.761900 0.001628 1.315352 18.260737 21.608761 8.465824];


colmat = cbrewer('qual','Paired',5);


if strcmp(type,'sips_sips')
    
    for i=1:numseries
        [out, term1, term2] = model_SipsSips(params(i).series,RHvec);
        plot(RHvec,out,'Color',colmat(i,:),'LineWidth',2)
        plot(RHvec,term1,'Color',colmat(i,:),'LineWidth',2)
        plot(RHvec,term2,'Color',colmat(i,:),'LineWidth',2)
    end
    
    [out, term1, term2] = model_SipsSips(opt_fit,RHvec);
    plot(RHvec,out,'Color',colmat(numseries + 1,:),'LineWidth',2)
    plot(RHvec,term1,'Color',colmat(numseries + 1,:),'LineWidth',2)
    plot(RHvec,term2,'Color',colmat(numseries + 1,:),'LineWidth',2)
        
end

if strcmp(type,'Do')
    
    for i=1:numseries
        plot(RHvec,model_DoLang(params(i).series,RHvec),'Color',colmat(i,:),'LineWidth',2)
    end
    
end


if isempty(leginfo) ~= true
    [rw,cl] = size(leginfo);
    
    if cl > rw 
        leginfo = leginfo';
    end
    
end

lcell = [{'exp #: 140828';'static rest';'static fit to 140828'}; leginfo];
legh = legend(lcell,'Location','North');
set(legh,'FontSize',10)
set(legh,'Box','off')
% print('-depsc','isotherms.eps')
export_fig 'isotherms.eps' 'isotherms.png'

function xi = model_DoLang(parw,x)

S0 = parw(1); k = parw(2);
Smu = parw(3); Kmu = parw(4); m = parw(5);
xi = S0*k*x./(1 + k*x)...
    + Smu*Kmu*x.^m./(1 + Kmu*x.^m);
end

function [xi,term1,term2] = model_SipsSips(parw,x)

S0 = parw(1); k = parw(2); m0 = parw(3);
Smu = parw(4); Kmu = parw(5); m = parw(6);

term1 = S0*k*x.^m0./(1 + k*x.^m0);
term2 = Smu*Kmu*x.^m./(1 + Kmu*x.^m);

xi = term1 + term2;

end



end