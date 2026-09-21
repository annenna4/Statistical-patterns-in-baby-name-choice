% Simulation of both versions of the age-structured neutral model described
% in the manuscript 'Innovation, conservatism, and the need to be distinguishable shape the diversity of first names' by
% Anne Kandler, Rafael D'Andrea, James O'Dwyer

clear
close all

nPop = 10^6; % population size
pMut = 3*10^-2; % innovation rate (per transmission event)
pDeath = 0.015; % per capita death rate
lambda = 0.4; % strength of anti-novelty bias

copyAll = 0; % 0: age-constrained copy pool as determined by copyThresholdHigh and copyThresholdLow defined below
% 1: copying from the whole population
copyThresholdHigh = 2; % upper bound of the age of the copying pool (ATTENTION: coincides with c_thresh+1 in the manuscript)
copyThresholdLow = 0; % lower bound of the age of the copying pool

tMax = 100; % time steps to be run after equilibrium has been reached
itMax = 10; % number of simulations

localMode = 0; % 0: age-structured transmission model as explained in section 2.1.,
% 1: age-structured transmission model with local interactions as explained in section 2.2
groupSize = 100; % size of local groups when localMode =1

PDmode = 0; % 0: do nothing, 1: calculate progeny distribution
saveMode = 0; % 0: do nothing, 1: save simulation output 

if ~isfolder('data')
    mkdir('data');
    fprintf('Created output folder: data\n');
end

% pre-allocation of space
coeffPL = zeros(itMax,4);
coeffPLpop = zeros(itMax,3);
mu_est = zeros(1,itMax);
if saveMode == 1
    cAll = cell(1,itMax);
    cPDAll = cell(1,itMax);
    lifeTimeAll = cell(1,itMax);
end

for sim = 1:itMax

    fprintf('Simulation # %d\n',sim)

    % burn-in period for reaching stationarity
    fprintf('Burn-in period\n')
    lifetimes = zeros(1,10^4);
    [pop,age] = get_burnIn(pDeath,nPop,pMut,copyAll,copyThresholdHigh,copyThresholdLow,lambda,localMode,groupSize,lifetimes);

    % re-name variant types for convenience
    names = unique(pop(1,:)); % unique names in pop
    [~, pop(1,:)] = ismember(pop(1,:),names);
    [~, idx_age] = ismember(names,age(1,:));
    h = age(2,idx_age);
    value = numel(names);
    valueIni = value;
    tini = max(pop(2,:));

    % pre-allocation of space
    maxVariants = 10^7; % guess of the maximum number of variants generated throughout the simulation
    age = zeros(2,maxVariants);
    age_count = numel(names);
    age(1,1:age_count) = 1:age_count;
    age(2,1:age_count) = h;
    death = zeros(1,numel(h));
    namesFreq = zeros(1,maxVariants);

    % calculating the cultural composition of the population for tMax time
    % steps
    fprintf('Generating populations \n')

    for t = tini+1:tini+tMax
        t-tini
        if localMode == 0 % age-structured neutral model as explained in section 2.1.
            [pop,value,namesFreq,age,lifetimes,age_count] = get_dynamics(t,pop,value,pDeath,nPop,pMut,lambda,copyAll,copyThresholdHigh,copyThresholdLow,PDmode,namesFreq,age,lifetimes,age_count);
        elseif localMode == 1 % age-structured neutral model with local interactions as explained in section 2.2
            [pop,value,namesFreq,addMut,age,lifetimes,age_count] = get_dynamics_local(t,pop,value,pDeath,nPop,copyAll,copyThresholdHigh,copyThresholdLow,PDmode,namesFreq,groupSize,lambda,age,lifetimes,age_count);
            effMut(t-tini) = addMut; % innovation rate in time step t
        end
    end

    % calculation of innovation rate of local model version
    if localMode == 1
        mu_est(sim) = mean(effMut);
    end

    fprintf('Calculating VAD \n')
    [types, ~, ic] = unique(pop(1,:)); % variant types present in population
    freq = accumarray(ic, 1)'; % frequencies of those types
    fprintf('Estimation of powerlaw - VAD \n')
    [alpha,xmin,L] = plfit(freq,'range',[1.01:0.01:3.01]); % using estimator developed by Clauset et al. (2009)
    coeffPLpop(sim,:) = [alpha,xmin,L];

    % plotting VAD
    [types, ~, ic] = unique(freq);
    c = accumarray(ic, 1) ./ length(freq);
    c = [[types'; types(end)+1], 1-[0; cumsum(c)]];
    c(c(:,2) < 10^-10, :) = [];
    maxFreqPop(sim) = c(end,1)./nPop; % frequency of most common variant type

    figure(1)
    loglog(c(:,1),c(:,2),'r--'); hold on;
    title('Variant abundance distribution')
    xlabel('Abundance x')
    ylabel('Probability P(x) of number of variants with abundance >= x')

    if saveMode == 1
        lifeTimeAll{sim} = lifetimes;
        cAll{sim} = c;
        if copyAll == 0
            name = sprintf('data/VAD_N%02d_pDeath%02d_cThresh%02d_lambda%02d.txt',nPop,pDeath,copyThresholdHigh,lambda);
            save(name,'-ASCII','coeffPLpop');
            name = sprintf('data/VADFULL_N%02d_pDeath%02d_cThresh%02d_lambda%02d.mat',nPop,pDeath,copyThresholdHigh,lambda);
            save(name,'cAll')
            name = sprintf('data/lifetime_N%02d_pDeath%02d_cThresh%02d_lambda%02d.mat',nPop,pDeath,copyThresholdHigh,lambda);
            save(name,'lifeTimeAll')
            if localMode == 1
                name = sprintf('data/effMut__N%02d_pDeath%02d_cThresh%02d_lambda%02d.txt',nPop,pDeath,copyThresholdHigh,lambda);
                save(name,'-ASCII','mu_est');
            end
        elseif copyAll == 1
            name = sprintf('data/VAD_N%02d_pDeath%02d_ALL_lambda%02d.txt',nPop,pDeath,lambda);
            save(name,'-ASCII','coeffPLpop');
            name = sprintf('data/VADFULL_N%02d_pDeath%02d_ALL_lambda%02d.mat',nPop,pDeath,lambda);
            save(name,'cAll')
            name = sprintf('data/lifetime_N%02d_pDeath%02d_ALL_lambda%02d.mat',nPop,pDeath,lambda);
            save(name,'lifeTimeAll')
            if localMode == 1
                name = sprintf('data/effMut__N%02d_pDeath%02d_ALL_lambda%02d.txt',nPop,pDeath,lambda);
                save(name,'-ASCII','mu_est');
            end
        end
    end

    if PDmode == 1 % calculation of progeny distribution

        fprintf('Calculating PD\n')
        namesFreq = nonzeros(namesFreq);
        namesFreq = reshape(namesFreq,numel(namesFreq),1);
        sumProgeny = sum(namesFreq);

        fprintf('Estimation of powerlaw - PD \n')
        [alpha,xmin,L] = plfit(namesFreq,'range',[1.01:0.01:3.01]); % using estimator developed by Clauset et al. (2009)
        coeffPL(sim,1:3) = [alpha,xmin,L];
        coeffPL(sim,4) = sumProgeny;

        % plotting PD
        [types, ~, ic] = unique(namesFreq);
        c1 = accumarray(ic, 1);
        c = c1./sum(c1);
        c = [[types; types(end)+1], 1-[0; cumsum(c)]];
        c(c(:,2) < 10^-10, :) = [];

        figure(2)
        loglog(c(1:end,1),c(1:end,2),'b','LineWidth',1); hold on;
        title('Progeny distribution')
        xlabel('Abundance x')
        ylabel('Probability P(x) of number of variants with abundance >= x')

        if saveMode == 1
            cPDAll{sim} = c;
            if copyAll == 0
                name = sprintf('data/PDcoeff_N%02d_pDeath%02d_cThresh%02d_lambda%02d_set%02d.txt',nPop,pDeath,copyThresholdHigh,lambda);
                save(name,'-ASCII','coeffPL');
                name = sprintf('data/PDFULL_N%02d_pDeath%02d_cThresh%02d_lambda%02d_set%02d.mat',nPop,pDeath,copyThresholdHigh,lambda);
                save(name,'cPDAll')
            elseif copyAll == 1
                name = sprintf('check/PDcoeff_N%02d_pDeath%02d_ALL_lambda%02d_set%02d.txt',nPop,pDeath,lambda);
                save(name,'-ASCII','coeffPL');
                name = sprintf('check/PDFULL_N%02d_pDeath%02d_ALL_lambda%02d_set%02d.mat',nPop,pDeath,lambda);
                save(name,'cPDAll')
            end
        end
    end
end

