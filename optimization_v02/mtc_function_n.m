%% mass transfer scaling 
% contains a number of (empirical) models specified in'type' for the mass 
% transfer scaling function in terms of the relative humidity
% outputs a <figure_name_type>.png figure where 'type' is the model

function mtc = mtc_function_n(mtc0,mtcscale,typee,parw)

rv = linspace(0,1,1e3);
nofrv = @(rv) model_SipsSips(parw,rv)/model_SipsSips(parw,1);

if nargin < 3
	typee = 'custom';
end

switch typee
    
    case 'monotonic'
        
        scaler = @(rr)  1 - mtcscale(1)*nofrv(rr).^mtcscale(2)./...
            (1 + mtcscale(1)*nofrv(rr).^mtcscale(2));        
    
    case 'monoplus'

        scaler = @(rr)  1 + mtcscale(1)*mtcscale(2)*nofrv(rr).^mtcscale(3)./...
            (1 + mtcscale(2)*nofrv(rr).^mtcscale(3));
    
    case 'exponential'

        scaler = @(rr) mtcscale(1)*exp(-mtcscale(2)*nofrv(rr)) + 1;
    
    case 'custom'
        
        scaler = @(rr) mtcscale(1)*exp(-mtcscale(2)*nofrv(rr))+ mtcscale(3)*mtcscale(4)*nofrv(rr).^mtcscale(5)./...
            (1 + mtcscale(4)*nofrv(rr).^mtcscale(5)) + 1;
end
           

mtc = mtc0*scaler(rv);

if nargin < 1
		   
	figure(1)
	clf
	hold on
	plot(nofrv(rv)*model_SipsSips(parw,1),scaler(rv),'k')
	xlabel('n [mol/kg]')
	ylabel('scaling factor')
else
	figure(1)
	clf
	hold on
	plot(nofrv(rv)*model_SipsSips(parw,1),mtc,'k')
	xlabel('n [mol/kg]')
	ylabel('k_{H_2O} [1/s]')
end
	
	
export_fig(strcat('mtc_func_',typee,'.png'))


figure(1) 
clf 
hold on 
set(gca,'yscale','log')
plot(nofrv(rv)*model_SipsSips(parw,1),mtc,'k')
xlabel('n [mol/kg]')
ylabel('mtc [1/s]')

export_fig(strcat('mtc_func_',typee,'_logy.png'))




function [xi,term1,term2] = model_SipsSips(parw,x)

S0 = parw(1); k = parw(2); m0 = parw(3);
Smu = parw(4); Kmu = parw(5); m = parw(6);

term1 = S0*k*x.^m0./(1 + k*x.^m0);
term2 = Smu*Kmu*x.^m./(1 + Kmu*x.^m);

xi = term1 + term2;

end

end





