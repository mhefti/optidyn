%% mass transfer scaling 
% contains a number of (empirical) models specified in'type' for the mass 
% transfer scaling function in terms of the relative humidity
% outputs a <figure_name_type>.png figure where 'type' is the model
rv = linspace(0,1,1e3);
type = 'monoplus';

switch type
    
    case 'monotonic'
        mtcscale(1) = [20];
        mtcscale(2) = [15];
        scaler = @(rr)  1 - mtcscale(1)*rr.^mtcscale(2)./...
            (1 + mtcscale(1)*rr.^mtcscale(2));        
    
    case 'monoplus'
        mtcscale(1) = [9];
        mtcscale(2) = [50000];
        mtcscale(3) = 15;
        mtcscale = [99;422.479300000000;17.5831400000000]
        scaler = @(rr)  1 + mtcscale(1)*mtcscale(2)*rr.^mtcscale(3)./...
            (1 + mtcscale(2)*rr.^mtcscale(3));
    
    case 'exponential'
        mtcscale(1) = [90];
        mtcscale(2) = [25.5];
        scaler = @(rr) mtcscale(1)*exp(-mtcscale(2)*rr) + 1;
    
    case 'custom'
        mtcscale(1) = [90];
        mtcscale(2) = 25.5;
        mtcscale(3) = 80;
        mtcscale(4) = 100;
        mtcscale(5) = 10;
        
        mtcscale6 = linspace(0,1,10);

        scaler = @(rr) mtcscale(1)*exp(-mtcscale(2)*rr) + mtcscale(3)*mtcscale(4)*rr.^mtcscale(5)./...
            (1 + mtcscale(4)*rr.^mtcscale(5)) + 1;
end
           

figure(1)
clf
hold on
plot(rv,scaler(rv),'k')
xlabel('rh')
ylabel('scaling factor')
export_fig(strcat('mtc_func_',type,'.png'))
