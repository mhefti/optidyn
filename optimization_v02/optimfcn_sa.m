function res = optimfcn_sa(param,y,hand,execname,expinfo,brut,initvals)
% optimfcn_sa is the objective function that is minimized by the
% 'algorithm' chosen in the input of optimization_dynamic_v0X; res is the
% value of the objective function, here the maximum likelihood estimate.
% The current parameters (absolute) and res are logged in fitting_log.txt
% input: 
%           param       vector of normalized parameters
%           hand        function handle pointing to the fitness function
%           y           structure containing the experimental data
%           execname    string w/o extension .exe
%           expinfo     structure containing information conc experiments
%           brut        logical flag to indicate usage of brutus
%           initvals    vector of initial guesses
% output: 
%           res         scalar w value of objective function


% param is normalized
param_abs = param.*initvals;

% get the model output 
md = hand(param_abs,execname,expinfo,brut,initvals);

% maximum likelihood
res_y_sum = zeros(expinfo.numexp,1);
res_T_sum = zeros(expinfo.numexp,1);

% constructs a string 'q_DATE_mode' where DATE is the date of the
% experiment specified in expinfo.list and mode is adsorption (a) or
% desorption (d) specified in expinfo.modes. 
% q is 'yH2O' or 'T05'
sfield = @(q,ii) sprintf('%s_%s_%s',q,num2str(expinfo.list(ii)),...
    expinfo.modes(ii));


% function handle for linear interpolation
% q is 'yH2O' or 'T05'
intercust = @(q,ii) interp1(md.(sfield('time',ii)),md.(sfield(q,ii)),...
    y.(sfield('time',ii)),'spline');


switch expinfo.phitype
    
    case 'MLE' % gives the maximum likelihood estimate
        
        for i = 1:expinfo.numexp
    
            if expinfo.fity
                ssq.(sfield('yH2O',i)) = sum( ( ...
                    y.(sfield('yH2O',i))/expinfo.maxy(i) - ...
                    intercust('yH2O',i)/expinfo.maxy(i) ...
                    ).^2);
                res_y_sum(i) = ssq.(sfield('yH2O',i));
            end
            
            if expinfo.fitT
                ssq.(sfield('T05',i)) = sum( ( ...
                    y.(sfield('T05',i))/expinfo.maxT(i) - ...
                    intercust('T05',i)/expinfo.maxT(i) ...
                    ).^2);
                res_T_sum(i) = ssq.(sfield('T05',i));
            end
            
        end
        
        
    case 'DV' % gives the derivative weighted MLE 
        
        for i = 1:expinfo.numexp
            
            % smoothing function handle
            sfunc = @(iny) smooth(iny,200); % change # to control degree of
                                            % smoothing
            
            % define generic function to compute the (positive) weights
            wfunc = @(iny,inx) abs((iny(2:end) - iny(1:end-1))./...
                (inx(2:end) - inx(1:end-1)));
            
            if expinfo.fity
                yv_mod = intercust('yH2O',i);
                yv_exp = sfunc(y.(sfield('yH2O',i)));
                tv_exp = y.(sfield('time',i));
                
                % compute the normalized weights
                w = wfunc(yv_exp,tv_exp);
                w = w/max(w);
                
                % skip some points to prevent weird effects at boundaries
                iskip = 20;
                eskip = 20;
                w = w(iskip:end-eskip);
                
                ssq.(sfield('yH2O',i)) = sum( w'*( ...
                    yv_exp(iskip:end-eskip-1)/expinfo.maxy(i) - ...
                    yv_mod(iskip:end-eskip-1)/expinfo.maxy(i) ...
                    ).^2);
                res_y_sum(i) = ssq.(sfield('yH2O',i));
            end
            
            if expinfo.fitT
                Tv_mod = intercust('T05',i);
                Tv_exp = sfunc(y.(sfield('T05',i)));
                tv_exp = y.(sfield('time',i));
                
                % compute the normalized weights
                w = wfunc(Tv_exp,tv_exp);
                w = w/max(w);
                
                % skip some points to prevent weird effects at boundaries
                iskip = 20;
                eskip = 20;
                w = w(iskip:end-eskip);
                
                ssq.(sfield('T05',i)) = sum( w'*( ...
                    Tv_exp(iskip:end-eskip-1)/expinfo.maxT(i) - ...
                    Tv_mod(iskip:end-eskip-1)/expinfo.maxT(i) ...
                    ).^2);
                res_T_sum(i) = ssq.(sfield('T05',i));
            end
            
        end

end


fprintf('-------------------fitting info---------------------\n')

for i = 1:expinfo.numexp
fprintf('error in experiment                    %s\n',num2str(expinfo.list(i)))

    if expinfo.fity
        fprintf('composition:          log(ssq(yH2O))   %f\n',log(res_y_sum(i)))
    end
    
    if expinfo.fitT
        fprintf('temperature:          log(ssq(T05))    %f\n',log(res_T_sum(i)))
    end
    
end
    
% if both fity and fitT are enabled, the log-likelihood is: 
% res = sum_{i=1}^(Nexp) log(ssq_y_i) + sum_{i=1}^(Nexp) log(ssq_T_i)
% in this way every experiment is equally weighted
if expinfo.fity && expinfo.fitT
    res = sum(log(res_y_sum)) + sum(log(res_T_sum));    
elseif expinfo.fity
    res = sum(log(res_y_sum));
else
    res = sum(log(res_T_sum));    
end

% write log file 
if expinfo.fitiso % if fitting isotherm enabled, write also isotherm pars
    customwrite('fitting_log.txt',{[param_abs' md.iso_pars res]},'\t','a+')
else
    customwrite('fitting_log.txt',{[param_abs' res]},'\t','a+')
end

fprintf('current value of objective function    %f\n',res)


