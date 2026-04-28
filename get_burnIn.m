
function [pop1,age1] = get_burnIn(pDeath,nPop,pMut,copyAll,copyThresholdHigh,copyThresholdLow,lambda,localMode,binSize,lifetimes)

% BURN IN RULE: it is assumed that the stationary state has been reached when all initially present variant types have gone extinct, ie. when all variant types have undergone neutral
% dynamic. This is a relatively strict rule that potentially takes a large amount of time
% 
% ageIni = ceil(rand(1,nPop)*50); % initial birth years
% t = max(ageIni);
% 
% maxVariants = 10^7;
% 
% % initialisation of population
% numbTypeIni = 2; % number of types initially present
% value = numbTypeIni;
% pop1 = get_popIni(nPop,numbTypeIni,ageIni);
% age1 = zeros(2,maxVariants); % innovation times of variants 
% age_count = numbTypeIni;
% age1(1,1:age_count) = 1:age_count;
% age1(2,1:age_count) = ones(1,age_count);
% 
% while  min(pop1(1,:))<numbTypeIni+1 % until all types have undergone neutral dynamics
%     t = t+1;
%     if localMode == 0
%         [pop1,value,~,age1,~,age_count] = get_dynamics(t,pop1,value,pDeath,nPop,pMut,lambda,copyAll,copyThresholdHigh,copyThresholdLow,0,[],age1,lifetimes,age_count);
%     elseif localMode == 1
%        [pop1,value,~,~,age1,~,age_count] = get_dynamics_local(t,pop1,value,pDeath,nPop,copyAll,copyThresholdHigh,copyThresholdLow,0,[],binSize,lambda,age1,lifetimes,age_count);
%     end
% end

% ALTERNATIVE BURN IN RULE: here it is assumed that that stationary state has been reached if two the heterogeneity index of two populations, both initialised with opposing initial conditions
% (one with maximum heterogeneity, one with minimum heterogeneity), have crossed.

ageIni = ceil(rand(1,nPop)*50); % initial birth years
t = max(ageIni);
tini = t;

maxVariants = 10^7;

initialisation of population 1
numbTypeIni = 2; % number of types initially present
value1 = numbTypeIni;
pop1 = get_popIni(nPop,numbTypeIni,ageIni);
age1 = zeros(2,maxVariants); % innovation times of variants 
age_count1 = numbTypeIni;
age1(1,1:age_count1) = 1:age_count1;
age1(2,1:age_count1) = ones(1,age_count1);

initialisation of population 2
numbTypeIni = nPop/10; % number of types initially present
value2 = numbTypeIni;
pop2 = get_popIni(nPop,numbTypeIni,ageIni);
age2 = zeros(2,maxVariants); % innovation times of variants 
age_count2 = numbTypeIni;
age2(1,1:age_count2) = 1:age_count2;
age2(2,1:age_count2) = ones(1,age_count2);

diffPop1Pop2 = 5;
stepTime = 1000;

while diffPop1Pop2>0
    for t = tini:tini
    t = t+1;

    if localMode == 0
        [pop1,value1,~,age1,~,age_count1] = get_dynamics(t,pop1,value1,pDeath,nPop,pMut,lambda,copyAll,copyThresholdHigh,copyThresholdLow,0,[],age1,lifetimes,age_count1);
        [pop2,value2,~,age2,~,age_count2] = get_dynamics(t,pop2,value2,pDeath,nPop,pMut,lambda,copyAll,copyThresholdHigh,copyThresholdLow,0,[],age2,lifetimes,age_count2);
    elseif localMode == 1
        [pop1,value1,~,~,age1,~,age_count1] = get_dynamics_local(t,pop1,value1,pDeath,nPop,copyAll,copyThresholdHigh,copyThresholdLow,0,[],binSize,lambda,age1,lifetimes,age_count1);
        [pop2,value2,~,~,age2,~,age_count2] = get_dynamics_local(t,pop2,value2,pDeath,nPop,copyAll,copyThresholdHigh,copyThresholdLow,0,[],binSize,lambda,age2,lifetimes,age_count2);
    end                                                       

    if mod(t,stepTime) == 0
        type = unique(pop1(1,:));
        h = hist(pop1(1,:),type)./nPop;
        divPop1 = sum(h.^2);
        type = unique(pop2(1,:));
        h = hist(pop2(1,:),type)./nPop;
        divPop2 = sum(h.^2);
        diffPop1Pop2 = divPop1-divPop2;
        if diffPop1Pop2 > 0.1
            stepTime = 500;
        elseif diffPop1Pop2 <0.1 && diffPop1Pop2 > 0.0
            stepTime = 100;
        end
    end
end
