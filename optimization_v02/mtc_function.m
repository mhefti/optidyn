%% mass transfer scaling 
% contains a number of (empirical) models specified in'type' for the mass 
% transfer scaling function in terms of the relative humidity
% outputs a <figure_name_type>.png figure where 'type' is the model

function mtc = mtc_function(mtc0,mtcscale,typee)

rv = linspace(0,1,1e3);

if nargin < 3
	typee = 'custom';
end

switch typee
    
    case 'monotonic'
        
        scaler = @(rr)  1 - mtcscale(1)*rr.^mtcscale(2)./...
            (1 + mtcscale(1)*rr.^mtcscale(2));        
    
    case 'monoplus'

        scaler = @(rr)  1 + mtcscale(1)*mtcscale(2)*rr.^mtcscale(3)./...
            (1 + mtcscale(2)*rr.^mtcscale(3));
    
    case 'exponential'

        scaler = @(rr) mtcscale(1)*exp(-mtcscale(2)*rr) + 1;
    
    case 'custom'
        
        scaler = @(rr) mtcscale(1)*exp(-mtcscale(2)*rr) + mtcscale(3)*mtcscale(4)*rr.^mtcscale(5)./...
            (1 + mtcscale(4)*rr.^mtcscale(5)) + 1;
end
           


mtc = mtc0*scaler(rv);

% if nargin < 1
% 		   
% 	figure(1)
% 	clf
% 	hold on
% 	plot(rv,scaler(rv),'k')
% 	xlabel('rh')
% 	ylabel('scaling factor')
% else
% 	figure(1)
% 	clf
% 	hold on
% 	plot(rv,mtc,'k')
% 	xlabel('rh')
% 	ylabel('mtc [1/s]')
% end
% % 	
% 	
% export_fig(strcat('mtc_func_',typee,'.png'))
% 
% 
% figure(1) 
% clf 
% hold on 
% set(gca,'yscale','log')
% plot(rv,mtc,'k')
% xlabel('rh')
% ylabel('mtc [1/s]')
% 
% export_fig(strcat('mtc_func_',typee,'_logy.png'))
% 



