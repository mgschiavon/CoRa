%% Figure S3. .
%   Gómez-Schiavon & El-Samad
%   August 2020
clear

%% Panel A
sim.an = 'ExSSs';
sim.ex = 'FigS2E';
sim.pp = 'mU';          % Label for perturbation type
sim.ax = 'mU';          % Label for condition/range
myXL = 'U synthesis rate (\mu_U)';
myYL = 'CoRa_{\mu_U\in\Theta}(\mu_U)';

fig = figure();
fig.Units = 'inches';
fig.PaperPosition = [2 1 4 3]*0.75;
fig.Position = fig.PaperPosition;
hold on;
    load(cat(2,'DATA_',sim.an,'_ATFv1_',sim.ex,'_',sim.pp,'_',sim.ax,'.mat'))
    plot(rho.values,DYs,'DisplayName','ATF v1',...
        'LineWidth',3,'Color',[1 0.6 0.78])
    load(cat(2,'DATA_',sim.an,'_ATFv2_',sim.ex,'_',sim.pp,'_',sim.ax,'.mat'))
    plot(rho.values,DYs,'DisplayName','ATF v2',...
        'LineWidth',3,'LineStyle','--','Color',[0.5 0.18 0.56])
        ylabel(myYL)
        ylim([0 1])
        xlabel(myXL)
        xlim([min(rho.values) max(rho.values)])
        box on
        set(gca,'XScale','log',...
            'XTick',[0.0001 0.01 1 100],'XGrid','on','XMinorGrid','off')
        legend
        
        print(gcf,cat(2,'RAW_',sim.ex,'_',sim.pp,'_',sim.ax),'-dpng','-r300')
        
sim.an = 'ExSSs';
sim.ex = 'FigS2E';
sim.pp = 'mU';          % Label for perturbation type
sim.ax = 'mY';          % Label for condition/range
myXL = 'Y synthesis rate (\mu_Y)';
myYL = 'CoRa_{\mu_Y\in\Theta}(\mu_U)';

fig = figure();
fig.Units = 'inches';
fig.PaperPosition = [2 1 4 3]*0.75;
fig.Position = fig.PaperPosition;
hold on;
    load(cat(2,'DATA_',sim.an,'_ATFv1_',sim.ex,'_',sim.pp,'_',sim.ax,'.mat'))
    plot(rho.values,DYs,'DisplayName','ATF v1',...
        'LineWidth',3,'Color',[1 0.6 0.78])
    load(cat(2,'DATA_',sim.an,'_ATFv2_',sim.ex,'_',sim.pp,'_',sim.ax,'.mat'))
    plot(rho.values,DYs,'DisplayName','ATF v2',...
        'LineWidth',3,'LineStyle','--','Color',[0.5 0.18 0.56])
        ylabel(myYL)
        ylim([0 1])
        xlabel(myXL)
        xlim([min(rho.values) max(rho.values)])
        box on
        set(gca,'XScale','log',...
            'XTick',[0.0001 0.01 1 100],'XGrid','on','XMinorGrid','off')
        legend
        
        print(gcf,cat(2,'RAW_',sim.ex,'_',sim.pp,'_',sim.ax),'-dpng','-r300')
        

sim.an = 'ExSSs';
sim.ex = 'FigS2E';
sim.pp = 'mW';          % Label for perturbation type
sim.ax = 'mY';          % Label for condition/range
myXL = 'Y synthesis rate (\mu_Y)';
myYL = 'CoRa_{\mu_Y\in\Theta}(\mu_W)';

fig = figure();
fig.Units = 'inches';
fig.PaperPosition = [2 1 4 3]*0.75;
fig.Position = fig.PaperPosition;
hold on;
    load(cat(2,'DATA_',sim.an,'_ATFv1_',sim.ex,'_',sim.pp,'_',sim.ax,'.mat'))
    plot(rho.values,DYs,'DisplayName','ATF v1',...
        'LineWidth',3,'Color',[1 0.6 0.78])
    load(cat(2,'DATA_',sim.an,'_ATFv2_',sim.ex,'_',sim.pp,'_',sim.ax,'.mat'))
    plot(rho.values,DYs,'DisplayName','ATF v2',...
        'LineWidth',3,'LineStyle','--','Color',[0.5 0.18 0.56])
        ylabel(myYL)
        ylim([0 1])
        xlabel(myXL)
        xlim([min(rho.values) max(rho.values)])
        box on
        set(gca,'XScale','log',...
            'XTick',[0.0001 0.01 1 100],'XGrid','on','XMinorGrid','off')
        legend
        
        print(gcf,cat(2,'RAW_',sim.ex,'_',sim.pp,'_',sim.ax),'-dpng','-r300')

%% END