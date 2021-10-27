%%                CoRa ANALYSIS                 %%
%%% DATA: Import output files from  DY_Main.jl %%%
%%%       to Matlab to generate figures.       %%%
% Mariana Gómez-Schiavon
% July, 2019

clear;
sim.mm = 'ATFv1';       % Label for motif file
sim.ex = 'Ex01';        % Label for parameters file
sim.pp = 'mY';          % Label for perturbation type
sim.ax = 'mY';          % Label for condition/range
sim.an = 'DYms';       % Chose analysis type (Options: ExSSs, ExDyn, DYms, OptDY)

%% Load & parse data
x = importdata(cat(2,'OUT_',sim.an,'_',sim.mm,'_',sim.ex,'_',sim.pp,'_',sim.ax,'.txt'),'\t',1);
if(strcmp(sim.an,'ExSSs'))
    rho.cond   = x.textdata{1};
    rho.values = x.data(:,1);
    for i = 2:(size(x.data,2)-1)
        str = x.textdata{i};
            [Is,Ie] = regexp(str,"_.+");
        ss.(str(1:3)).(str((Is+1):Ie)) = x.data(:,i);
    end
    rho.pert = x.textdata{end};
    DYs = x.data(:,end);
    clear Is Ie i x str
elseif(strcmp(sim.an,'ExDyn'))
    ll = x.textdata;
    pF = x.data([x.data(:,1)==1],2);
    tF = x.data([x.data(:,1)==1],3);
    sF = x.data([x.data(:,1)==1],4:end);
    pN = x.data([x.data(:,1)==0],2);
    tN = x.data([x.data(:,1)==0],3);
    sN = x.data([x.data(:,1)==0],4:end);
    lN = x.textdata;
    clear x
elseif(strcmp(sim.an,'DYms'))
    for i = 1:size(x.data,2)
        str = x.textdata{i};
        [Is,Ie] = regexp(str,"\d+\.\d+");
        if(Is & Ie)
            break;
        else
            p.(str) = x.data(:,i);
        end
    end
    DYs = x.data(:,i:end);
    rho.name   = sim.pp;
    rhoV = x.textdata(1,i:end);
    rho.values = zeros(size(rhoV));
    for i = 1:length(rhoV)
        rho.values(i) = str2num(rhoV{i});
    end
    clear i rhoV Ie Is x str
elseif(strcmp(sim.an,'OptDY'))
    rho.name = sim.pp;
    if(size(x.data,2)>601)
        DYs.curv = x.data(:,(end-600):end);
        x.data = x.data(:,1:(end-601));
        x.textdata = x.textdata(:,1:(end-601));
    end
    epsT = x.textdata{end-1};
        [Is,Ie] = regexp(epsT,"\d\.\d+");
        epsT  = epsT(Is:Ie);
        clear Is Ie 
    R = x.data(:,1);
    if(length(R)==2*max(R))
        I = x.data(:,2);
        DYi.min = x.data([I == 0],end);
        DYi.thr = x.data([I == 0],end-1);
        DYf.min = x.data([I ~= 0],end);
        DYf.thr = x.data([I ~= 0],end-1);
        for i = 3:(size(x.data,2)-2)
            pi.(x.textdata{i}) = x.data([I == 0],i);
            pf.(x.textdata{i}) = x.data([I ~= 0],i);
        end
    else
        R0 = 1;
        while(R0<=max(R))
            I = x.data(:,2);
            DYs(R0).min = x.data([R == R0],end);
            DYs(R0).thr = x.data([R == R0],end-1);
            for i = 3:(size(x.data,2)-2)
                ps(R0).(x.textdata{i}) = x.data([R == R0],i);
            end
            R0 = R0 + 1;
        end
    end
    clear i x R I R0
else
    'ERROR: Undetermined analysis. Options: ExSSs, DYms, OptD'
end
save(cat(2,'DATA_',sim.an,'_',sim.mm,'_',sim.ex,'_',sim.pp,'_',sim.ax,'.mat'))

%% FIGURE
if(strcmp(sim.an,'ExSSs'))
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 3 3];
    fig.Position = fig.PaperPosition;

    plot(rho.values,DYs,'DisplayName','CoRa(\mu_Y)',...
        'LineWidth',3,'Color',[0 0 0])
        xlabel('Y synthesis rate (\mu_Y)','FontSize',12)
        xlim([0.001 1000])
        ylabel('CoRa_{\mu_Y\in\Theta}(\mu_Y)','FontSize',12)
        ylim([0 1])
        set(gca,'XScale','log','XTick',10.^[-2:2:2],'FontSize',12)
        box on    
elseif(strcmp(sim.an,'ExDyn'))
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
        title(cat(2,'CoRa_{\mu_Y\in\Theta}(\mu_Y)=',...
            num2str(log10(sF(end,1)/sF(1,1))/log10(sN(end,1)/sN(1,1)),2)))
        set(gca,'YTick','','YMinorTick','Off')
        set(gca,'YScale','log',...
            'YTick',(10.^[0:(log10(1.05)/3):log10(1.05)])*sF(1,1),... 'YTickLabel',{cat(2,'logY|\mu_Y^{(',sim(i).No,')}'),'','','','',''},...        'YTickLabelRotation',80,...
            'XTick',[500:8000:10000],'XTickLabel',[500:8000:10000]-500,...
            'XGrid','on','YGrid','on')
        box on
elseif(strcmp(sim.an,'DYms'))
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 3 3];
    fig.Position = fig.PaperPosition;
    C = [0 1 1;0 0.962745070457458 1;0 0.925490200519562 1;0 0.888235330581665 1;0 0.850980401039124 1;0 0.813725471496582 1;0 0.776470601558685 1;0 0.739215731620789 1;0 0.701960802078247 1;0 0.664705872535706 1;0 0.627451002597809 1;0 0.596078455448151 1;0 0.564705908298492 1;0 0.533333361148834 1;0 0.501960813999176 1;0 0.470588266849518 1;0 0.439215689897537 1;0 0.407843142747879 1;0 0.376470595598221 1;0 0.345098048448563 1;0 0.313725501298904 1;0 0.285205006599426 0.909090936183929;0 0.256684511899948 0.818181812763214;0 0.228164002299309 0.727272748947144;0 0.199643507599831 0.636363625526428;0 0.171122997999191 0.545454561710358;0 0.142602503299713 0.454545468091965;0 0.114082001149654 0.363636374473572;0 0.0855614989995956 0.272727280855179;0 0.0570410005748272 0.181818187236786;0 0.0285205002874136 0.0909090936183929;0 0 0;0 0 0;0.0909090936183929 0.0142602501437068 0.0909090936183929;0.181818187236786 0.0285205002874136 0.181818187236786;0.272727280855179 0.0427807494997978 0.272727280855179;0.363636374473572 0.0570410005748272 0.363636374473572;0.454545468091965 0.0713012516498566 0.454545468091965;0.545454561710358 0.0855614989995956 0.545454561710358;0.636363625526428 0.0998217537999153 0.636363625526428;0.727272748947144 0.114082001149654 0.727272748947144;0.818181812763214 0.128342255949974 0.818181812763214;0.909090936183929 0.142602503299713 0.909090936183929;1 0.156862750649452 1;1 0.191721141338348 1;1 0.226579532027245 1;1 0.261437922716141 1;1 0.296296298503876 1;1 0.331154674291611 1;1 0.366013079881668 1;1 0.400871455669403 1;1 0.43572986125946 1;1 0.470588237047195 1;1 0.499108731746674 1;1 0.527629256248474 1;1 0.55614972114563 1;1 0.58467024564743 1;1 0.613190710544586 1;1 0.641711235046387 1;1 0.670231759548187 1;1 0.698752224445343 1;1 0.727272748947144 1;1 0.755793213844299 1;1 0.7843137383461 1];
    hold on;
        pp = fieldnames(p); pp = pp{1};
        pN = regexprep(pp,'^m','\\mu_'); pN = regexprep(pN,'^g','\\gamma_'); pN = regexprep(pN,'^e','\\eta_'); pN = regexprep(pN,'_P','_+'); pN = regexprep(pN,'_M','_-');pN = regexprep(pN,'^b','\\beta_');pN = regexprep(pN,'kD','K_D');pN = regexprep(pN,'_$','');pN = regexprep(pN,'_Us','_\{Us\}');
        rhoN = regexprep(rho.name,'^m','\\mu_'); rhoN = regexprep(rhoN,'^g','\\gamma_'); rhoN = regexprep(rhoN,'^e','\\eta_'); rhoN = regexprep(rhoN,'_P','_+'); rhoN = regexprep(rhoN,'_M','_-');rhoN = regexprep(rhoN,'^b','\\beta_');rhoN = regexprep(rhoN,'kD','K_D');rhoN = regexprep(rhoN,'_$','');rhoN = regexprep(rhoN,'_Us','_\{Us\}');
        for j = 1:7
            plot(rho.values,DYs(j,:),'DisplayName','CoRa(\mu_Y)',...
                'LineWidth',3,'Color',C(1+((j-1)*10),:))
        end
            xlabel('Y synthesis rate (\mu_Y)','FontSize',12)
            xlim([0.001 1000])
            ylabel('CoRa_{\mu_Y\in\Theta}(\mu_Y)','FontSize',12)
            ylim([0 1])
            if(strcmp(pp,'mW'))
                title('Changing W synthesis rate','FontSize',12)
            elseif(strcmp(pp,'mU'))
                title('Changing U synthesis rate','FontSize',12)
            end
            set(gca,'XScale','log','XTick',10.^[-2:2:2],'FontSize',12)
            box on

            annotation('textbox',[0.6 0.3 0.25 0.125],...
                'String',{cat(2,pN,'= ',num2str(p.(pp)(4)))},...
                'FitBoxToText','on','LineWidth',3,'FontSize',12);
end
        print(gcf,cat(2,'RAW_',sim.an,'_',sim.mm,'_',sim.ex,'_',sim.pp,'_',sim.ax),'-dpng','-r300')
%% END