%% Toy model: cumulative normal nonlinearity with reasonably rectifying parameterization
% center & surround "activations" (i.e. post-filter generator signals)
% combined in either an independent or shared -type cascade
% Measure rectification index as in data
clear all; close all; clc;

xx = (-10:0.1:10);
alpha = 1;
beta = 0.15;
gamma = -1.1;
epsilon = -0.15;

NL=normcdfNonlinearity(xx,alpha,beta,gamma,epsilon);

figure(3); clf; plot(xx, NL,'k-')

%center, surround activations
cc = linspace(-5,5,200);
ss = linspace(-5,5,200);
[CC, SS] = meshgrid(cc,ss);
[~, zeroInd] = min(abs(cc)); % find "Zero" response index


%indep model. NL then sum
resp_indep = normcdfNonlinearity(CC(:),alpha,beta,gamma,epsilon) +...
    normcdfNonlinearity(SS(:),alpha,beta,gamma,epsilon);
%shared NL model. Sum then NL
resp_shared = normcdfNonlinearity(CC(:) + SS(:),alpha,beta,gamma,epsilon); 

resp_indep = reshape(resp_indep,length(cc),length(ss));
resp_shared = reshape(resp_shared,length(cc),length(ss));

colors = pmkmp(length(cc));
figure(1); clf;
subplot(221)
surf(cc,ss,resp_indep)
subplot(222); hold on;
RI_indep = [];
for ll = 1:length(cc)
    tempResp = resp_indep(:,ll);
    rPlus = tempResp(end);
    rMinus = tempResp(1);
    r0 = tempResp(zeroInd);
    RI_indep(ll) = 1 - (r0 - rMinus)./(rPlus - r0);
    plot(cc,tempResp,'Color',colors(ll,:)) 
end

subplot(223)
surf(cc,ss,resp_shared)
subplot(224); hold on;
RI_shared = [];
for ll = 1:length(cc)
    tempResp = resp_shared(:,ll);
    rPlus = tempResp(end);
    rMinus = tempResp(1);
    r0 = tempResp(zeroInd);
    RI_shared(ll) = 1 - (r0 - rMinus)./(rPlus - r0);
    
    plot(cc,tempResp,'Color',colors(ll,:)) 
    
end

figure; clf; fig1=gca; initFig(fig1,'Surround activation','Rect. index')
addLineToAxis(ss,RI_indep,'indep',fig1,[0.5 0.5 0.5],'-','none')
addLineToAxis(ss,RI_shared,'shared',fig1,'k','-','none')
makeAxisStruct(fig1,'rectIndexToyModel' ,'RFSurroundFigs')


%% Same thing for piecewise linear, rectifying nonlinearity
clear all; close all; clc;

%center, surround activations
cc = -10:1:10;
ss = -10:1:10;
[CC, SS] = meshgrid(cc,ss);

resp_indep = max(CC(:),0) + max(SS(:),0); %indep model. NL then sum
resp_shared = max(CC(:) + SS(:),0); %shared NL model. Sum then NL

resp_indep = reshape(resp_indep,length(cc),length(ss));
resp_shared = reshape(resp_shared,length(cc),length(ss));

colors = pmkmp(length(cc));
figure(1); clf;
subplot(221)
surf(cc,ss,resp_indep)
subplot(222); hold on;
RI_indep = [];
for ll = 1:length(cc)
    tempResp = resp_indep(:,ll);
    rPlus = tempResp(end);
    rMinus = tempResp(1);
    r0 = tempResp(11);
    RI_indep(ll) = 1 - (r0 - rMinus)./(rPlus - r0);
    plot(cc,tempResp,'Color',colors(ll,:)) 
end

subplot(223)
surf(cc,ss,resp_shared)
subplot(224); hold on;
RI_shared = [];
for ll = 1:length(cc)
    tempResp = resp_shared(:,ll);
    rPlus = tempResp(end);
    rMinus = tempResp(1);
    r0 = tempResp(11);
    RI_shared(ll) = 1 - (r0 - rMinus)./(rPlus - r0);
    
    plot(cc,tempResp,'Color',colors(ll,:)) 
    
end

figure(2); clf; hold on;
plot(ss,RI_indep,'b')
plot(ss,RI_shared,'k')

