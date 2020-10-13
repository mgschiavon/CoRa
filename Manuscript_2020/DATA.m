%%                CoRa ANALYSIS                 %%
%%% DATA: Import output files from  DY_Main.jl %%%
%%%       to Matlab to generate figures.       %%%
% Mariana Gómez-Schiavon
% July, 2019

clear;
sim.mm = 'ATFv1';       % Label for motif file
sim.ex = 'Fig3';        % Label for parameters file
sim.pp = 'mY';          % Label for perturbation type
sim.ax = 'mY';          % Label for condition/range
sim.an = 'DYms';       % Chose analysis type (Options: ExSSs, DYms, OptDY)

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
delete(cat(2,'OUT_',sim.an,'_',sim.mm,'_',sim.ex,'_',sim.pp,'_',sim.ax,'.txt'))

%% END