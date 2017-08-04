% Claire Guerin - 03/08/2017
% Human Evolution Project in Montpellier
% Human Evolutionary Biology Team - ISEM - Bernard Godelle & Michel Raymond

% This script analyses the adaptation to nutrition data generated under
% Python

PathName = 'C:\Users\Claire\Dropbox\MEME\Montpellier\AdaptiveDynamicsStratification\pData';
cd(PathName)

fileList = dir([PathName,'\*.mat']);
nGen = 200000;
nStrat = 100;

nFiles = size(fileList);

fileName = fileList(nFiles).name;
lastGens = load(fileName);
testEqui = squeeze(lastGens.prob(end,:,:));

%%

XX = linspace(0,1,size(testEqui,2));
YY = linspace(0,1,size(testEqui,1));

fig1=figure(1);clf;
set(fig1,'defaulttextinterpreter','latex','Color','w');

contourf(X,Y,testEqui)
colormap(summer)
axis equal
xlabel('Resident diet strategy $\hat{x}$','fontsize',20)
ylabel('Mutant diet strategy $x$','fontsize',20)
caxis([0,1])
hcb = colorbar;
set(hcb,'defaulttextinterpreter','latex');
title(hcb,'p','fontsize',18)
%%
export_fig '..\dietEqui.bmp' -m2

any(any(testEqui == 1)) % = 0: there is never a situation where the population is only composed of residents.

%%

[X, Y] = find(testEqui>0);
nPolyM = numel(X);
subset = round(linspace(1,nPolyM,100));
percPolyM = nPolyM/numel(testEqui);

fig2 = figure(2);clf;
set(fig2,'defaulttextinterpreter','latex','Color','w');
xlim([0,2000])
ylim([0,1])
hold on

colorshades = summer(nPolyM);

threshLow = [0,.25,.75];
threshHigh = [.25,.75,1];

graphID = nan([1,nPolyM]);

for file = 1:1
    
    fileName = fileList(file).name;
    Stmp = load(fileName);
    generation = strsplit(fileName(5:end-4),'-');
    start = str2double(generation{1});
    nd = str2double(generation{2});
    
%     evolution = nan([1,nPolyM]);
    
    for i = 1:nPolyM
        
        resident = X(i);
        mutant = Y(i);
        
        evolution = squeeze(Stmp.prob(:,resident,mutant));
        meanP = nanmean(evolution);
        graphID(i) = find(meanP >= threshLow & meanP < threshHigh);
        
%         subplot(numel(threshLow),1,graphID)
        plot(start+1:nd,evolution,'-','col',colorshades(i,:))
        clear evolution mutant resident
%         evolution(i) = nanmean(squeeze(Stmp.prob(:,resident,mutant)));  
    end
    
%     histogram(evolution);
    clear Stmp fileName generation start nd
end

hold off