%% diffusion barrier - scaling function
% contains a number of (empirical) models specified in'type' for the  
% diffusion barrier
% outputs a <figure_name_type>.png figure where 'type' is the model
rv = linspace(0,1,1e3);
type = 'monoder';

switch type
    
    case 'monoder'
        
        s(1) = [50];
        s(2) = [1.1];
        s(3) = 1;
        % define the function
        % note: this function is derived as follows:
        % d/dx (k*x^m/(1 + k*x^m)) = (k*m*x^(-m - 1))/(1 + k*x^m)^2
        % modify the function empirically: 
        % ff(x) = (k*m*x^m)/(1 + k*x^m)^2
        % d/dx ff(x) = 0 -> solve for x_max: x_max = k^(-1/m)
        % finally: 
        % der = ff(x)/ff(x_max) -> bounded in [0,1]
        der = @(rr) 4*(s(1)*rr.^s(2))./(1 + s(1)*rr.^s(2)).^2;
        scaler = @(rr) 1 - s(3)*der(rr) ;
        
    case 'normal'
        
        s(1) = 0.1;
        s(2) = 2;
        s(3) = 30;
        s(4) = 0.8;
        % define the function
        phi = @(rr) 1/(2*pi)*exp(-s(3)*(rr - s(1)).^s(2)/2);
        scaler = @(rr) s(4)*phi(rr)/max(phi(rr));        
    
    case 'skewed_normal'
        
        s(1) = 0;
        s(2) = 2;
        s(3) = 1;
        % define the function
        phi = @(rr) 1/(2*pi)*exp(-(rr - s(1)).^2/2);
        bphi = @(rr) 1/2*(1 + erf((rr - s(1))/sqrt(2)));
        scaler = @(rr) s(3)*phi(rr).*bphi(s(2)*rr);
    
    case 'exponential'
        
        s(1) = 1;
        s(2) = 1;
        s(3) = 0.01;
        % define the function
        exfunc = @(rr) s(1)*(1 - rr).^s(2).*(exp(rr.^s(3)) - 1);
        scaler = @(rr) exfunc(rr);
        
    case 'quadratic'
        
        s(1) = 1;
        s(2) = 4;
        s(3) = 0.4;
        % define the function
        qfunc = @(rr) s(1)*(s(3) - rr).^s(2) + 1;
        scaler = @(rr) qfunc(rr);
            
    case 'tangens'
        
        s(1) = 1;
        s(2) = 20;
        s(3) = 0;
        % define the function
        tfunc = @(rr) s(1)*(0 + tanh(s(2)*rr + s(3)));
        scaler = @(rr) tfunc(rr);
                   
end

figure(1)
clf
hold on
plot(rv,scaler(rv),'k')
xlabel('rh')
ylabel('scaling factor')
export_fig(strcat('diff_barrier_',type,'.png'))