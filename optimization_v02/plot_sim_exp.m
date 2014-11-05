clear all
% close all


% use correlation for axial conductivity


%--- set which simulations to show along which experiments ----------------
% system(['.' filesep 'Fixed_Bed'])

%specify wether to recalculate all the specified simulations, just plot the
%already present data, or to recalculate just specific ones (specify
%further down):  'all', 'some' , 'none'
recalc = 'none'; 
plot_mode = 2;  %0 - no plot, 1 - plot, 2 - plot + pdf
% print2pdf = 0;

%--- Figure 1 ------------------------------------
 fig_name{1} = 'GS_36_ip';

%      [   sim  exp  ]
show{1}= [  36	79 ];
%             2	21  ;
%             3	20  ;
%             4	23  ;
%             5	27  ;
%             6	31  ];
        
        
% % % --- Figure 1 ------------------------------------
%  fig_name{1} = '10CO2_90N2_25C';
% 
% %      [   sim  exp  ]
% show{1}= [  1	30  ;
%             2	21  ;
%             3	20  ;
%             4	23  ;
%             5	27  ;
%             6	31  ];

% % --- Figure 2 ------------------------------------
% fig_name{2} = '10CO2_90N2_45C';
% % 
% %        [ sim  exp  ]
% show{2}= [  7	71  ;
%             8	69  ;
%             9	66  ;
%             10	67  ;
%             11	68  ;
%             12	64  ];
%         
% % --- Figure 3 ------------------------------------
% fig_name{3} = '50CO2_50N2_25C';
% % % 
% % %        [ sim  exp  ]
% show{3}= [  13	34  ;
%             14	14  ;
%             15	11   ;
%             16	77  ;
%             17	56  ;
%             18	75  ];
% %         
% % 
% %--- Figure 4 ------------------------------------
% fig_name{4} = '50CO2_50N2_45C';
% % 
% %        [ sim  exp  ]
% show{4}= [  19	73  ;
%             20	48  ;
%             21	40  ;
%             22	72  ;
%             23	44  ;
%             24	37  ];
% 
% %--- Figure 5 ------------------------------------
% fig_name{5} = '80CO2_20N2_25C';
% % 
% %        [ sim  exp  ]
% show{5}= [  25	92  ;
%             26	90  ;
%             27	83  ;
%             28	84  ;
%             29	85  ;
%             30	78  ];        
% % 
% % % %--- Figure 6 ------------------------------------
% fig_name{6} = '80CO2_20N2_45C';
% % % 
% % %        [ sim  exp  ]
% show{6}= [  31	89  ;
%             32	88  ;
%             33	81  ;
%             34	86  ;
%             35	87  ;
%             36	79  ];


%--------------------------------------------------------------------------
%   plotting parameters
%--------------------------------------------------------------------------
sw=44;  % screen width in cm
lw=1;   % line width
fs=12;  % font size
Msize = 3; % marker size

% define symbols for experimental points
sym = ['o';'o';'o';'o';'o';'o';'o';'o';'o';'o'];

% define colors for plots
t1 =   [0.0 0.6 0.0];  %  [0.8 0.2 0.2];  % colors for T10 - T110
t2 =   [0.5 0.9 0.0];  %  [0.8 0.7 0.2];  %
t3 =   [1.0 0.8 0.2];  %  [0.2 0.6 0.2];  %
t4 =   [1.0 0.0 0.0];  %  [0.2 0.2 0.9];  %
t5 =   [0.5 0.2 0.0];  %  [0.6 0.2 0.6];  %
c_n2 = [0.0 0.0 0.0];  %  [0.2 0.2 0.2];  % color for each component
c_co2= [0.0 0.0 1.0];  %  [0.3 0.3 0.8];  %
c_he = [0.0 1.0 0.0];  %  [0.2 0.7 0.2];  %

%% -----------------------  ALL FIGURES  ----------------------------------
Nfig = length(fig_name); %number of figures to be created

calc = [];

if (strcmp(recalc,'all'))
    for n=1:Nfig
        calc = [calc;show{n}(:,1)];
    end
elseif (strcmp(recalc,'some'))
    calc = [3]'; % define which experiments are recalculated here
%     calc = [7:12]'; % define which experiments are recalculated here
%     calc = [7 8]'; % define which experiments are recalculated here
%     calc = [28:36]'; % define which experiments are recalculated here
end

for n=1:length(calc)
    fid = fopen('iterations.dat','r+');
    fprintf(fid, '%u %u   // first/last experiment to evaluate',calc(n),calc(n));
    fclose(fid);
    system(['.' filesep 'Fixed_Bed']);
end

if (plot_mode >0)
for m = 1:Nfig;
    %--- collect all the necessary filepaths for figure m -----------------
    % Simulations
    Nsim = length(show{m}(:,1));
    simfile_Y = cell(Nsim,1);
    simfile_Y_pipe = cell(Nsim,1);
    simfile_T = cell(Nsim,1);

    for i=1:Nsim
        if(show{m}(i,1)<10)
            simfile_Y{i}     =['.\experiment_0' int2str(show{m}(i,1)) '\exitprofile.txt'];
            simfile_Y_pipe{i}=['.\experiment_0' int2str(show{m}(i,1)) '\exitprofile.txt'];
            simfile_T{i}     =['.\experiment_0' int2str(show{m}(i,1)) '\temperatures.txt'];
        else
            simfile_Y{i}     =['.\experiment_' int2str(show{m}(i,1)) '\exitprofile.txt'];
            simfile_Y_pipe{i}=['.\experiment_' int2str(show{m}(i,1)) '\exitprofile.txt'];
            simfile_T{i}     =['.\experiment_' int2str(show{m}(i,1)) '\temperatures.txt'];
        end
    end
    
    % Experiments
    N_comp = 3;
    fid = fopen('T:\adsorption\column\Measurements\AP3-60_N2\exp_list.csv');
    exp_list = textscan(fid, '%d %f %f %s %s %*[^\n]','delimiter',',','Headerlines',1);
    fclose(fid);
    exp_T = exp_list{1,2};
    exp_P = exp_list{1,3};
    LV_files = exp_list{1,4};
    MS_files = exp_list{1,5};
    show_exp = show{m}(:,2);
    LV = LV_files(show_exp);
    MS = MS_files(show_exp);
    
    %--- read the data---------------------------------------------------------
    % Simulations
    mf_sim_prof=cell(Nsim,1); mf_sim_pipe_prof=cell(Nsim,1); 
    T_sim_prof=cell(Nsim,1); T_sim_C=cell(Nsim,1);
    for j=1:Nsim
        mf_sim_prof{j} = dlmread(simfile_Y{j},'',1,0);
        mf_sim_pipe_prof{j} = dlmread(simfile_Y_pipe{j},'',1,0);
        T_sim_prof{j} = dlmread(simfile_T{j},'',1,0);
        T_sim_C{j} = [T_sim_prof{j}(:,1) (T_sim_prof{j}(:,2:end)-273.15)];
    end

    % Experiments
    [res] = read_all(LV,MS,N_comp);
    % res is a cell array with a row for each experiment. Each row contains:
    % -----------------------------------------------------------------------
    % ¦ time (LV) ¦ flow ¦  P  ¦  T  ¦ time(MS) ¦  y  ¦ MS signal ¦ species ¦
    % ¦   (L,1)   ¦ (L,1)¦(L,3)¦(L,5)¦   (M,1)  ¦(M,n)¦   (M,n)   ¦  (n,1)  ¦
    % -----------------------------------------------------------------------
    % where L is the number of measurement points done by LabView, M is the
    % number of points measured by the MS, and n is the number of components.
    % Therefore, the call res{2,4}(150,3) returns the 150th temperature
    % measurement for TI 3 for experiment 2
    %--------------------------------------------------------------------------
    
    %--- make the figure ------------------------------------------------------

    figure('Name',fig_name{m},'NumberTitle','off')
    if(Nsim>4)
        fw=sw;
    else
        fw=2+(Nsim*(sw-2)/5);
    end
    set(gcf,'units','centimeters','Position',[0.2 3 fw 22])
    
    minT=zeros(1,Nsim); maxT=zeros(1,Nsim);

    
    for k=1:Nsim;
            
        Tinit = mean(res{k,4}(1,:),2);
        avgP  = mean(mean(res{k,3}(:,1:2),2),1);
        avgF  = mean(res{k,2});
        cond = {['T= ' int2str(Tinit) ' °C'];
                ['P= ' int2str(avgP) ' bar'];
                ['Flow= ' int2str(avgF) ' cm^3/s']};

        [res_mf] = residuals(mf_sim_pipe_prof{k},[res{k,5} res{k,6}],T_sim_C{k},[res{k,1} res{k,4}], 'mle');
        resid = {'residuals:';
               ['exit= ' num2str(res_mf)]};
 
        subplot('Position',[(2/fw+((1-(2/fw))/Nsim)*(k-1)) 0.6 ((1-(2/fw))/Nsim)*0.9 0.3])
        if(k==1)
            set(gca,'Fontsize',fs,'ytick',[0:0.2:1])
            ylabel('mole fractions [-]', 'Fontsize',fs)
        else
            set(gca,'Fontsize',fs,'ytick',[0:0.2:1],'YTickLabel',{'','','','','',''})
        end
        box on
        grid on
        xlimits=[0 max(mf_sim_pipe_prof{k}(end,1), res{k,5}(end,1))];
        xlim(xlimits)
        ylim([0 1.02])
        xlabel('time [s]', 'Fontsize',fs)
        hold on
%         text(0.7*xlimits(2),0.9,'N_2','FontSize',fs,'Color',c_n2)
%         text(0.7*xlimits(2),0.1,'CO_2','FontSize',fs,'Color',c_co2)
        text(0,-0.3,cond,'Fontsize',fs)
        plot(mf_sim_pipe_prof{k}(:,1),mf_sim_pipe_prof{k}(:,2),'-','Color',c_co2,'LineWidth',lw)
        plot(mf_sim_pipe_prof{k}(:,1),mf_sim_pipe_prof{k}(:,3),'-','Color',c_n2,'LineWidth',lw)
        % plot(mf_sim_prof{k}(:,1),mf_sim_prof{k}(:,4),'-','Color',c_he,'LineWidth',lw)
        plot(res{k,5},res{k,6}(:,1),sym(1),'Color',c_n2,'MarkerSize',Msize)
        plot(res{k,5},res{k,6}(:,2),sym(1),'Color',c_co2,'MarkerSize',Msize)
        % plot(res{k,5},res{k,6}(:,3),sym(1),'Color',c_he,'MarkerSize',Msize)
        hold off
        

%         minT(k) = ;
%         maxT(k) = max(max(res{k,4}));

        mintemp=10*floor(0.1*min(min(res{k,4})));
        maxtemp=10*ceil(0.1*max(max(res{k,4})));
        subplot('Position',[(2/fw+((1-(2/fw))/Nsim)*(k-1)) 0.15 ((1-(2/fw))/Nsim)*0.88 0.3])
        if(k==1)
            set(gca,'Fontsize',fs,'ytick',[0:10:100])
            ylabel('Temperature [°C]', 'Fontsize',fs)
        else
            set(gca,'Fontsize',fs,'ytick',[0:10:100])
        end
%         set(gca,'Fontsize',fs,'ytick',[0:10:100])
        box on
        grid on
        xlim([0 max(T_sim_C{k}(end,1), res{k,1}(end,1))])
        ylim([mintemp maxtemp])
        xlabel('time [s]', 'Fontsize',fs)
%         ylabel('Temperature [°C]', 'Fontsize',fs)
        hold on
        plot(res{k,1},res{k,4}(:,1),sym(2),'Color',t1,'MarkerSize',Msize)
        plot(res{k,1},res{k,4}(:,2),sym(2),'Color',t2,'MarkerSize',Msize)
        plot(res{k,1},res{k,4}(:,3),sym(2),'Color',t3,'MarkerSize',Msize)
        plot(res{k,1},res{k,4}(:,4),sym(2),'Color',t4,'MarkerSize',Msize)
        plot(res{k,1},res{k,4}(:,5),sym(2),'Color',t5,'MarkerSize',Msize)
        plot(T_sim_C{k}(:,1),T_sim_C{k}(:,2),'-','Color',t1,'LineWidth',lw)
        plot(T_sim_C{k}(:,1),T_sim_C{k}(:,3),'-','Color',t2,'LineWidth',lw)
        plot(T_sim_C{k}(:,1),T_sim_C{k}(:,4),'-','Color',t3,'LineWidth',lw)
        plot(T_sim_C{k}(:,1),T_sim_C{k}(:,5),'-','Color',t4,'LineWidth',lw)
        plot(T_sim_C{k}(:,1),T_sim_C{k}(:,6),'-','Color',t5,'LineWidth',lw)

        hold off
        [maxT, imax] = max(res{k,4}(:,1:5));
%         [maxT, imax] = max(T_sim_C{k}(:,2:6));
%         text(res{k,1}(imax(1)),maxtemp-4,'10','FontSize',fs,'Color',t1)
%         text(res{k,1}(imax(2)),maxtemp-4,'35','FontSize',fs,'Color',t2)
%         text(res{k,1}(imax(3)),maxtemp-4,'60','FontSize',fs,'Color',t3)
%         text(res{k,1}(imax(4)),maxtemp-4,'85','FontSize',fs,'Color',t4)
%         text(res{k,1}(imax(5)),maxtemp-4,'110 cm','FontSize',fs,'Color',t5)
%         text(T_sim_C{k}(imax(1),1),maxT(1)+3,'10','FontSize',fs,'Color',t1)
%         text(T_sim_C{k}(imax(2),1),maxT(2)+3,'35','FontSize',fs,'Color',t2)
%         text(T_sim_C{k}(imax(3),1),maxT(3)+3,'60','FontSize',fs,'Color',t3)
%         text(T_sim_C{k}(imax(4),1),maxT(4)+3,'85','FontSize',fs,'Color',t4)
%         text(T_sim_C{k}(imax(5),1),maxT(5)+3,'110 cm','FontSize',fs,'Color',t5)
%         text(0, mintemp-0.3*(maxtemp-mintemp),resid,'BackgroundColor',[1 1 1],'EdgeColor','k')
        
    end
    
    if(plot_mode == 2)        
        set(gcf,'PaperPositionMode','auto')
        print('-depsc', [fig_name{m} '.eps'])
        system(['epstopdf ' fig_name{m} '.eps']);
    end
    
end
end