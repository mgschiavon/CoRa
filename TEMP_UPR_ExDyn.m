sim.mm = 'UPR';       % Label for motif file
sim.ex = 'Ex03d3';        % Label for parameters file
sim.pp = 'cD';          % Label for perturbation type
sim.ax = 'cD';          % Label for condition/range
sim.an = 'ExDyn';       % Chose analysis type (Options: ExSSs, DYms, OptDY)

x = importdata(cat(2,'OUT_',sim.an,'_',sim.mm,'_',sim.ex,'_',sim.pp,'_',sim.ax,'.txt'),'\t',1);
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
    plot(tF,sF(:,1)+sF(:,3)+sF(:,4),'LineWidth',2,'Color',[1 0.6 0.78])
    plot(tN,sN(:,1)+sN(:,3)+sN(:,4),'LineWidth',2,'LineStyle','--','Color',[0 0 0])
    xlabel('Time (min)')
%     xlim(500+([-.1 1]*8000))
%     ylim((10.^[-(log10(1.05)/10) log10(1.05)+(log10(1.05)/10)])*sF(1,1))
    ylabel('Output')
    title(cat(2,'DTT=',num2str(3000/(0.001516231919278*6.022e+05*2.15434934),3),...
        'mM | CoRa_{c_D\cdotDTT\in\Theta}(c_D\cdotDTT)=',...
        num2str(log10(sum(sF(end,[1,3,4]))/sum(sF(1,[1,3,4])))/log10(sum(sN(end,[1,3,4]))/sum(sN(1,[1,3,4]))),2)))
    set(gca,'YTick','','YMinorTick','Off')
    set(gca,'YScale','log',...
        'YTick',(10.^[0:(log10(1.05)/3):log10(1.05)])*sF(1,1),... 'YTickLabel',{cat(2,'logY|\mu_Y^{(',sim.No,')}'),'','','','',''},...        'YTickLabelRotation',80,...        'XTick',[500:8000:10000],'XTickLabel',[500:8000:10000]-500,...
        'XGrid','on','YGrid','on')
    box on
    print(gcf,cat(2,'RAW_',sim.ex,'_',sim.mm),'-dpng','-r300')