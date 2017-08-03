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
fig1=figure(1);clf;
set(fig1,'defaulttextinterpreter','latex','Color','w');

contourf(testEqui)
colormap(summer)
axis equal
xlabel('Resident diet strategy $\hat{x}$','fontsize',20)
ylabel('Mutant diet strategy $x$','fontsize',20)
hcb = colorbar;
set(hcb,'defaulttextinterpreter','latex');
title(hcb,'p','fontsize',18)
export_fig '..\dietEqui.bmp' -m2

fullDynamic = nan(nGen); % evolution of p for 1 strategy

for file = 1:nFiles
    fileName = fileList(file).name;
    Stmp = load(fileName);
    pos = strstrip(fileName(5:end-4),'-');
    start = char2num(pos{1});
    nd = char2num(pos{2});
    fullMat(start:nd,:,:) = Stmp.prob;
    clear Stmp
end


