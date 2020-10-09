%% Figure 2.  Characterizing the antithetic feedback motif (ATF) using CoRa
%   Gómez-Schiavon & El-Samad
%   July 2020
clear

%% Panel B
sim.an = 'ExSSs';
sim.ex = 'Fig2B';
sim.pp = 'mY';          % Label for perturbation type
sim.ax = 'mY';          % Label for condition/range
myXL = 'Y synthesis rate (\mu_Y)';
myYL = 'CoRa_{\mu_Y\in\Theta}(\mu_Y)';

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
        
        plot([0.3863 0.3863],[0 1],'LineWidth',3,'Color',[0.8 0.8 0.8]-0.1,...
            'DisplayName','\mu_Y^{(1)}')
        plot([3.9 3.9],[0 1],'LineWidth',3,'Color',[0.8 0.8 0.8],...
            'DisplayName','\mu_Y^{(2)}')
        plot([125 125],[0 1],'LineWidth',3,'Color',[0.8 0.8 0.8]+0.1,...
            'DisplayName','\mu_Y^{(3)}')
        
        print(gcf,cat(2,'RAW_',sim.ex),'-dpng','-r300')

%% Panel C
sim.an = 'ExSSs';
sim.ex = 'Fig2C';
sim.pp = 'mY';          % Label for perturbation type
sim.ax = 'eM';          % Label for condition/range
myXL = 'C degradation rate (\eta_-)';
myYL = 'CoRa_{\eta_-\in\Theta}(\mu_Y)';

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
        
        print(gcf,cat(2,'RAW_',sim.ex),'-dpng','-r300')

%% Panel D
clear sim
sim(1).mm = 'ATFv1';       % Label for motif file
sim(1).an = 'ExDyn';
sim(1).pp = 'mY';          % Label for perturbation type
sim(1).ax = 'mY';          % Label for condition/range
sim(1).ex = 'Fig2DL';
sim(1).No = '1';

sim(2).mm = 'ATFv2';       % Label for motif file
sim(2).an = 'ExDyn';
sim(2).pp = 'mY';          % Label for perturbation type
sim(2).ax = 'mY';          % Label for condition/range
sim(2).ex = 'Fig2DL';
sim(2).No = '1';

sim(3).mm = 'ATFv1';       % Label for motif file
sim(3).an = 'ExDyn';
sim(3).pp = 'mY';          % Label for perturbation type
sim(3).ax = 'mY';          % Label for condition/range
sim(3).ex = 'Fig2DM';
sim(3).No = '2';

sim(4).mm = 'ATFv2';       % Label for motif file
sim(4).an = 'ExDyn';
sim(4).pp = 'mY';          % Label for perturbation type
sim(4).ax = 'mY';          % Label for condition/range
sim(4).ex = 'Fig2DM';
sim(4).No = '2';

sim(5).mm = 'ATFv1';       % Label for motif file
sim(5).an = 'ExDyn';
sim(5).pp = 'mY';          % Label for perturbation type
sim(5).ax = 'mY';          % Label for condition/range
sim(5).ex = 'Fig2DR';
sim(5).No = '3';

sim(6).mm = 'ATFv2';       % Label for motif file
sim(6).an = 'ExDyn';
sim(6).pp = 'mY';          % Label for perturbation type
sim(6).ax = 'mY';          % Label for condition/range
sim(6).ex = 'Fig2DR';
sim(6).No = '3';

for i = 1:length(sim)
    x = importdata(cat(2,'OUT_',sim(i).an,'_',sim(i).mm,'_',sim(i).ex,'_',sim(i).pp,'_',sim(i).ax,'.txt'),'\t',1);
    ll = x.textdata;
    pF = x.data([x.data(:,1)==1],2);
    tF = x.data([x.data(:,1)==1],3);
    sF = x.data([x.data(:,1)==1],4:end);
    pN = x.data([x.data(:,1)==0],2);
    tN = x.data([x.data(:,1)==0],3);
    sN = x.data([x.data(:,1)==0],4:end);
    lN = x.textdata;
    clear x

    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 3 1.5];
    fig.Position = fig.PaperPosition;
    hold on;
        plot(tF,sF(:,1),'LineWidth',2,'Color',[1 0.6 0.78])
        plot(tN,sN(:,1),'LineWidth',2,'LineStyle','--','Color',[0 0 0])
        xlabel('Time (min)')
        xlim(500+([-.1 1]*8000))
        ylim((10.^[-(log10(1.05)/10) log10(1.05)+(log10(1.05)/10)])*sF(1,1))
        ylabel('Output')
        title(cat(2,'CoRa_{\mu_Y^{(',sim(i).No,')}\in\Theta}(\mu_Y)=',...
            num2str(log10(sF(end,1)/sF(1,1))/log10(sN(end,1)/sN(1,1)),2)))
        set(gca,'YTick','','YMinorTick','Off')
        set(gca,'YScale','log',...
            'YTick',(10.^[0:(log10(1.05)/3):log10(1.05)])*sF(1,1),... 'YTickLabel',{cat(2,'logY|\mu_Y^{(',sim(i).No,')}'),'','','','',''},...        'YTickLabelRotation',80,...
            'XTick',[500:8000:10000],'XTickLabel',[500:8000:10000]-500,...
            'XGrid','on','YGrid','on')
        box on
        print(gcf,cat(2,'RAW_',sim(i).ex,'_',sim(i).mm),'-dpng','-r300')
end;