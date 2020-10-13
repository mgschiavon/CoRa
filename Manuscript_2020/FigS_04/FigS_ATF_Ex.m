%% Figure S3. Effect of dilution on the antithetic feedback (ATF) control.
%   Gómez-Schiavon & El-Samad
%   July 2020
clear

%% Panel A
simT.ex = 'FigS3g';
simT.mm = 'ATFv1';       % Label for motif file
simT.an = 'ExSSs';
simT.pp = 'mY';          % Label for perturbation type
simT.ax = 'mY';          % Label for condition/range
myXL = 'Y synthesis rate (\mu_Y)';
myYL = 'CoRa_{\mu_Y\in\Theta}(\mu_Y)';
gG = [1e-3,1e-5,1e-7];

% Hard limits
gY  = 1;         % Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
mU  = 0.125;     % U synthesis rate dependent of Y (nM/min)
gU  = 0;         % U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
mW  = 0.1;       % W constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
gW  = 0;         % W degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
e0  = 0.0001;    % U:W dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
eP  = 0.0375;    % U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
eM  = 0.5;       % U:W (C) mutual anhilation/co-degradation (e.g. [0.03,2] 1/nM 1/min)

for i = 1:3
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 1.5]*0.75;
    fig.Position = fig.PaperPosition;
    axes('Parent',fig,'Position',[0.2 0.35 0.75 0.6])
    hold on;
        load(cat(2,'DATA_',simT.an,'_',simT.mm,'_',simT.ex,num2str(i),'_',simT.pp,'_',simT.ax,'.mat'))
        plot(rho.values,DYs,'LineWidth',3,...
            'DisplayName',cat(2,'\gamma = ',num2str(gG(i),1),' min^{-1}'))
        ylabel(myYL)
        ylim([0 1])
        xlabel(myXL)
        xlim([min(rho.values) max(rho.values)])
        box on
        set(gca,'XScale','log',...
            'XTick',10.^[-6:3:9],'XGrid','on','XMinorGrid','off')
        legend
        
        print(gcf,cat(2,'RAW_FigS3_g',num2str(i),'_',simT.mm),'-dpng','-r300')
        
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 1.5]*0.75;
    fig.Position = fig.PaperPosition;
    axes('Parent',fig,'Position',[0.2 0.35 0.75 0.6])
    hold on;
        g = gG(i);
        load(cat(2,'DATA_',simT.an,'_',simT.mm,'_',simT.ex,num2str(i),'_',simT.pp,'_',simT.ax,'.mat'))
        plot(rho.values,(ss.FbR.Y(1))*rho.values/min(rho.values),...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,(ss.FbR.Y(end))*rho.values/max(rho.values),...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,ss.FbR.Y,'LineWidth',2)
            xlabel(myXL)
            xlim([min(rho.values) max(rho.values)])
            ylabel('Y_{ss}')
%             ylim([0.0001 10000])
            box on
            set(gca,'XScale','log','YScale','log','YTick',10.^[-12:4:12],...
                'XTick',10.^[-6:3:9],'XGrid','on','XMinorGrid','off')

            print(gcf,cat(2,'RAW_FigS3_g',num2str(i),'_',simT.mm,'_Y'),'-dpng','-r300')
            
            
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 1.5]*0.75;
    fig.Position = fig.PaperPosition;
    axes('Parent',fig,'Position',[0.2 0.35 0.75 0.6])
    hold on;
        g = gG(i);
        load(cat(2,'DATA_',simT.an,'_',simT.mm,'_',simT.ex,num2str(i),'_',simT.pp,'_',simT.ax,'.mat'))
        plot([min(rho.values) max(rho.values)],[0 0]+(mW/(g+gW)),...
            'LineWidth',2,'Color',[0.8 0.8 0.8])
        plot(rho.values,(mW/mU)*(gY+g)./rho.values,...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,ss.FbR.W,'LineWidth',2)
            xlabel(myXL)
            xlim([min(rho.values) max(rho.values)])
            ylabel('W_{ss}')
            ylim([min(ss.FbR.W)/102 max(ss.FbR.W)*100])
            box on
            set(gca,'XScale','log','YScale','log','YTick',10.^[-12:4:12],...
                'XTick',10.^[-6:3:9],'XGrid','on','XMinorGrid','off')

            print(gcf,cat(2,'RAW_FigS3_g',num2str(i),'_',simT.mm,'_limW'),'-dpng','-r300')
        
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 1.5]*0.75;
    fig.Position = fig.PaperPosition;
    axes('Parent',fig,'Position',[0.2 0.35 0.75 0.6])
    hold on;
        g = gG(i);
        load(cat(2,'DATA_',simT.an,'_',simT.mm,'_',simT.ex,num2str(i),'_',simT.pp,'_',simT.ax,'.mat'))
        plot([min(rho.values) max(rho.values)],[0 0]+(mW/(g+gW)),...
            'LineWidth',2,'Color',[0.8 0.8 0.8])
        plot([min(rho.values) max(rho.values)],[0 0]+(mW/(g+gW+eM)),...
            'LineWidth',2,'Color',[0.8 0.8 0.8])
        plot(rho.values,(mW/mU)*(gY+g)./rho.values,...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,ss.FbR.W + ss.FbR.C,'LineWidth',2)
            xlabel(myXL)
            xlim([min(rho.values) max(rho.values)])
            ylabel('W_T_{ss}=W_{ss}+C_{ss}')
            ylim([min(ss.FbR.W + ss.FbR.C)/102 max(ss.FbR.W + ss.FbR.C)*100])
            box on
            set(gca,'XScale','log','YScale','log','YTick',10.^[-12:4:12],...
                'XTick',10.^[-6:3:9],'XGrid','on','XMinorGrid','off')

            print(gcf,cat(2,'RAW_FigS3_g',num2str(i),'_',simT.mm,'_limWT'),'-dpng','-r300')
        
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 1.5]*0.75;
    fig.Position = fig.PaperPosition;
    axes('Parent',fig,'Position',[0.2 0.35 0.75 0.6])
    hold on;
        g = gG(i);
        load(cat(2,'DATA_',simT.an,'_',simT.mm,'_',simT.ex,num2str(i),'_',simT.pp,'_',simT.ax,'.mat'))
        plot(rho.values,(ss.FbR.U(1))*rho.values/min(rho.values),...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,(ss.FbR.U(end))*rho.values/max(rho.values),...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,ss.FbR.U,'LineWidth',2)
            xlabel(myXL)
            xlim([min(rho.values) max(rho.values)])
            ylabel('U_{ss}')
%             ylim([0.0001 10000])
            box on
            set(gca,'XScale','log','YScale','log','YTick',10.^[-12:4:12],...
                'XTick',10.^[-6:3:9],'XGrid','on','XMinorGrid','off')
            
            print(gcf,cat(2,'RAW_FigS3_g',num2str(i),'_',simT.mm,'_limUT'),'-dpng','-r300')
end

%% Panel B
simT.ex = 'FigS3g';
simT.mm = 'ATFv2';       % Label for motif file
simT.an = 'ExSSs';
simT.pp = 'mY';          % Label for perturbation type
simT.ax = 'mY';          % Label for condition/range
myXL = 'Y synthesis rate (\mu_Y)';
myYL = 'CoRa_{\mu_Y\in\Theta}(\mu_Y)';
gG = [1e-3,1e-5,1e-7];

% Hard limits
gY  = 1;         % Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
mU  = 0.125;     % U synthesis rate dependent of Y (nM/min)
gU  = 0;         % U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
mW  = 0.1;       % W constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
gW  = 0;         % W degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
e0  = 0.0001;    % U:W dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
eP  = 0.0375;    % U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
eM  = 0.5;       % U:W (C) mutual anhilation/co-degradation (e.g. [0.03,2] 1/nM 1/min)

for i = 1:3
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 1.5]*0.75;
    fig.Position = fig.PaperPosition;
    axes('Parent',fig,'Position',[0.2 0.35 0.75 0.6])
    hold on;
        load(cat(2,'DATA_',simT.an,'_',simT.mm,'_',simT.ex,num2str(i),'_',simT.pp,'_',simT.ax,'.mat'))
        plot(rho.values,DYs,'LineWidth',3,...
            'DisplayName',cat(2,'\gamma = ',num2str(gG(i),1),' min^{-1}'))
        ylabel(myYL)
        ylim([0 1])
        xlabel(myXL)
        xlim([min(rho.values) max(rho.values)])
        box on
        set(gca,'XScale','log',...
            'XTick',10.^[-6:3:9],'XGrid','on','XMinorGrid','off')
        legend
        
        print(gcf,cat(2,'RAW_FigS3_g',num2str(i),'_',simT.mm),'-dpng','-r300')
        
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 1.5]*0.75;
    fig.Position = fig.PaperPosition;
    axes('Parent',fig,'Position',[0.2 0.35 0.75 0.6])
    hold on;
        g = gG(i);
        load(cat(2,'DATA_',simT.an,'_',simT.mm,'_',simT.ex,num2str(i),'_',simT.pp,'_',simT.ax,'.mat'))
        plot(rho.values,(ss.FbR.Y(1))*rho.values/min(rho.values),...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,(ss.FbR.Y(end))*rho.values/max(rho.values),...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,ss.FbR.Y,'LineWidth',2)
            xlabel(myXL)
            xlim([min(rho.values) max(rho.values)])
            ylabel('Y_{ss}')
%             ylim([0.0001 10000])
            box on
            set(gca,'XScale','log','YScale','log','YTick',10.^[-12:4:12],...
                'XTick',10.^[-6:3:9],'XGrid','on','XMinorGrid','off')

            print(gcf,cat(2,'RAW_FigS3_g',num2str(i),'_',simT.mm,'_Y'),'-dpng','-r300')
            
            
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 1.5]*0.75;
    fig.Position = fig.PaperPosition;
    axes('Parent',fig,'Position',[0.2 0.35 0.75 0.6])
    hold on;
        g = gG(i);
        load(cat(2,'DATA_',simT.an,'_',simT.mm,'_',simT.ex,num2str(i),'_',simT.pp,'_',simT.ax,'.mat'))
        plot([min(rho.values) max(rho.values)],[0 0]+(mW/(g+gW)),...
            'LineWidth',2,'Color',[0.8 0.8 0.8])
        plot(rho.values,(mW/mU)*(gY+g)./rho.values,...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,(ss.FbR.W(end))./(rho.values/max(rho.values)),...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,ss.FbR.W,'LineWidth',2)
            xlabel(myXL)
            xlim([min(rho.values) max(rho.values)])
            ylabel('W_{ss}')
            ylim([min(ss.FbR.W)/2 max(ss.FbR.W)*2])
            box on
            set(gca,'XScale','log','YScale','log','YTick',10.^[-12:4:12],...
                'XTick',10.^[-6:3:9],'XGrid','on','XMinorGrid','off')

            print(gcf,cat(2,'RAW_FigS3_g',num2str(i),'_',simT.mm,'_limW'),'-dpng','-r300')
        
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 1.5]*0.75;
    fig.Position = fig.PaperPosition;
    axes('Parent',fig,'Position',[0.2 0.35 0.75 0.6])
    hold on;
        g = gG(i);
        load(cat(2,'DATA_',simT.an,'_',simT.mm,'_',simT.ex,num2str(i),'_',simT.pp,'_',simT.ax,'.mat'))
        plot([min(rho.values) max(rho.values)],[0 0]+(mW/(g+gW)),...
            'LineWidth',2,'Color',[0.8 0.8 0.8])
        plot([min(rho.values) max(rho.values)],[0 0]+(mW/(g+gW+eM)),...
            'LineWidth',2,'Color',[0.8 0.8 0.8])
        plot(rho.values,(mW/mU)*(gY+g)./rho.values,...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,ss.FbR.W + ss.FbR.C,'LineWidth',2)
            xlabel(myXL)
            xlim([min(rho.values) max(rho.values)])
            ylabel('W_T_{ss}=W_{ss}+C_{ss}')
            ylim([min(ss.FbR.W + ss.FbR.C)/102 max(ss.FbR.W + ss.FbR.C)*100])
            box on
            set(gca,'XScale','log','YScale','log','YTick',10.^[-12:4:12],...
                'XTick',10.^[-6:3:9],'XGrid','on','XMinorGrid','off')

            print(gcf,cat(2,'RAW_FigS3_g',num2str(i),'_',simT.mm,'_limWT'),'-dpng','-r300')
        
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 1.5]*0.75;
    fig.Position = fig.PaperPosition;
    axes('Parent',fig,'Position',[0.2 0.35 0.75 0.6])
    hold on;
        g = gG(i);
        load(cat(2,'DATA_',simT.an,'_',simT.mm,'_',simT.ex,num2str(i),'_',simT.pp,'_',simT.ax,'.mat'))
        plot(rho.values,(ss.FbR.U(1))*rho.values/min(rho.values),...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,(ss.FbR.U(end))*rho.values/max(rho.values),...
            'LineWidth',2,'LineStyle','--','Color',[0.8 0.8 0.8])
        plot(rho.values,ss.FbR.U,'LineWidth',2)
            xlabel(myXL)
            xlim([min(rho.values) max(rho.values)])
            ylabel('U_{ss}')
%             ylim([0.0001 10000])
            box on
            set(gca,'XScale','log','YScale','log','YTick',10.^[-12:4:12],...
                'XTick',10.^[-6:3:9],'XGrid','on','XMinorGrid','off')
            
            print(gcf,cat(2,'RAW_FigS3_g',num2str(i),'_',simT.mm,'_limUT'),'-dpng','-r300')
end

%% END