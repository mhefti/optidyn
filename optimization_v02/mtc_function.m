%% mtc-mono

% figure(1)
% rv = linspace(0,1,1e3);
% mtcscale(1) = [9255.37800000000];
% mtcscale(2) = [7.74613790000000];
% scaler = @(rr) 1 - mtcscale(1)*rr.^mtcscale(2)./(1 + mtcscale(1)*rr.^mtcscale(2));
% plot(rv,scaler(rv),'k')
% xlabel('rH [-]')
% ylabel('scaling factor {\it s} [-]')
% export_fig 'scaling_mono.eps' 'scaling_mono.png'
% 
% 
% 
% 
% scaler = @(rr) 1 - mtcscale(1)*rr.^mtcscale(2)./(1 + mtcscale(1)*rr.^mtcscale(2)) ... + 
%                mtcscale(3)*exp(-mtcscale(4)*(rr - mtcscale(5)));



rv = linspace(0,1,1e3);
type = 'monoplus';
% type = 'custom';
% type = 'exponential';
% scaler = @(rr) 1 + mtcscale(1)*rr.^mtcscale(2)./(1 + mtcscale(1)*rr.^mtcscale(2));
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
% for i = 1:length(mtcscale6)
%     mtcscale(6) = mtcscale6(i);
%     scaler = @(rr) 1 + mtcscale(1)*exp(-mtcscale(2)*(rr - mtcscale(6))) + mtcscale(1) + ...
%         mtcscale(5)*mtcscale(3)*rr.^mtcscale(4)./(1 + mtcscale(3)*rr.^mtcscale(4));   
%     plot(rv,scaler(rv),'k')
%     
% end

plot(rv,scaler(rv),'k')
xlabel('rh')
ylabel('scaling factor')
export_fig(strcat('mtc_func_',type,'.png'))
% figure(1)
% rv = linspace(0,1,1e3);
% mtcscale(1) = [9255.37800000000];
% mtcscale(2) = [7.74613790000000];
% % exponential
% mtcscale(3) = 1;
% mtcscale(4) = 29;
% mtcscale(5) = 0.1;
% % second term
% % mtcscale(3) = 20;
% % mtcscale(4) = 6;
% % mtcscale(5) = 0.4;
% 
% scaler = @(rr)  - mtcscale(1)*rr.^mtcscale(2)./(1 + mtcscale(1)*rr.^mtcscale(2)) + ... ... + 
%                mtcscale(3)*exp(mtcscale(4)*(rr - mtcscale(5)));
% scaler = @(rr) 1 - mtcscale(1)*rr.^mtcscale(2)./(1 + mtcscale(1)*rr.^mtcscale(2)) + ... ... + 
%                  mtcscale(3)*tanh(mtcscale(4)*(rr - mtcscale(5))) + mtcscale(3);           
% % scaler = @(rr) 1 - mtcscale(1)*rr.^mtcscale(2)./(1 + mtcscale(1)*rr.^mtcscale(2)) + ... ... + 
% %                mtcscale(3)*(rr - mtcscale(5)).^mtcscale(4)./...
% %                (1 + mtcscale(3)*(rr - mtcscale(5)).^mtcscale(4)); 
% plot(rv,scaler(rv),'k')
% xlabel('rH [-]')
% ylabel('scaling factor {\it s} [-]')
% export_fig 'scaling_mono.eps' 'scaling_mono.png'           
%            
           
           
%  %% arcus tangens type
% % figure(1)
% % clf
% % a = 2.9;
% % b = [1/2.5 0.75 6];
% % c = [1.6 7.5 7.5];
% % y = linspace(0,0.1)*100;
% % mf = @(y,a,b,c) 3-b(1)*atan(b(2)*(y - b(2))) + c(1)*atan(c(2)*(y - c(3)));
% % plot(y,mf(y,a,b,c)/mf(0,a,b,c),'k','LineWidth',2)
% 
% 
% figure(
% 
% 
% 
% figure(1)
% clf
% a = 2.9;
% b = [1/2.5 0.75 6];
% c = [1.6 7.5 7.5];
% y = linspace(0,0.1)*100;
% isopars = [1459.29990000000 0.00122150280000000 0.890484450000000 ...
%             25 5.25624610000000 5.98209800000000];
% 
% addpath('isotherm_fit');        
% addpar.iso_type = 'Sips_Sips';
% addpar.isoflex = false;
% 
% n = @(rh) isotherm_func(isopars,rh,addpar);        
% nsat = n(1);
% rhv = linspace(0,1);
% sat = @(rh) n(rh)./nsat;
% mf = @(rh,a,b,c) cos(1.7*atan(sat(4*rh.^1.5)));
% 
% a = 1; b = 50; c = 1;
% mf = @(rh) a./sqrt(a + 10000*sat(rh).^c);
% facv = linspace(0.1,10,50);
% sf = linspace(0,1,300);
% facv=5;
% hold on
% % for i =1:length(facv)
% % plot(rhv,(1 - sat(rhv).^facv(i)),'k','LineWidth',2)
% % end
% nn = 10; 
% k = 10;
% kv = linspace(1,100,20);
% ex = linspace(0.1,10,20);
% hold on
% for i = 1:length(kv)
%     for j = 1:length(ex)
%         k = kv(i);
%         exx = ex(j);
%         plot(rhv,(1 - k*rhv.^exx./(1 + k*rhv.^exx)))
%     end
% end
% clf
% k(1) = 0.5;
% k(2) = 1000;
% exx(1) = 0.5;
% exx(2) = 10;
% w(1) = 1;
% w(2) = 1 - w(1);
% 
% kv1 = linspace(0.01,100,25);
% kv2 = kv1; 
% rr = 1:length(kv1);
% hold on
% wv = linspace(0,1,100);
% 
% for i = rr
%     for j = rr
%         exx(1) = kv1(i);
%         exx(2) = kv2(j);
%         plot(rhv,(1 - kv1(i)*rhv.^exx(1)./(1 + kv1(i)*rhv.^exx(1)) + ...
%               kv2(j)*(1 - rhv).^exx(2)./(1 + kv2(j)*(1 - rhv).^exx(2))))
%     end
% end
% 
% 
% 
% % potential
% 
% figure(2) 
% 
% plot(rhv, -(10*rhv).^exx(1) + (5*rhv).^exx(2))
% 
% 
% 
% 
% % % ylim([0 1])          
% %           
% %           
% % figure(2) 
% % clf
% % a = 1;
% % hold on 
% % a = 0;
% % for i = 1:50
% %     a = a + 0.01;
% %     exfunc = @(rhv,a) (1 + exp(-a*rhv)).^(-1);
% %     plot(rhv,1 - erf(a*rhv))
% % end
% 
% % p1 = linspace(1,10,10);
% % p2 = linspace(1,10000,10);
% % p3 = linspace(1,20,10);
% % hold on
% % for i = 1:10
% %     a = p1(i);
% %     for j=1:10
% %         b = p2(j);
% %         for k=1:10
% %             c = p3(k);
% %             mf = @(rh) a./sqrt(a + b*sat(rh).^c);
% %             plot(rhv,mf(rhv))
% %         end
% %     end
% % end
% 
% %% exponential type
% % isopars = [1459.29990000000 0.00122150280000000 0.890484450000000 ...
% %             25 5.25624610000000 5.98209800000000];
% % n = @(rh) model_SipSips(isopars,rh);        
% % nsat = n(1);
% % sat = @(rh) n(rh)./nsat;
% % 
% % 
% % p1 = linspace(1,10,10);
% % p2 = linspace(0,1,10);
% % p3 = linspace(0,1,10);
% % 
% % pp = [4 0.6 0.8];
% % rhv = linspace(0,1);
% % 
% % figure(2) 
% % clf 
% % hold on
% % for i = 1:10
% %     p1t = p1(i);
% %     for j=1:10
% %         p2t = p2(j);
% %         for k=1:10
% %             p3t = p3(k);
% %             pp = [p1t p2t p3t];
% %             mf = @(rh,pp) 1 + pp(2)*sat(rh).^pp(3).*(exp(pp(1)*sat(rh)) - 1);
% %             plot(rhv,mf(rhv,pp))
% %         end
% %     end
% % end