function [pop,value,namesFreq,addMut,age,lifetimes,age_count] = get_dynamics_local(t,pop,value,pDeath,nPop,copyAll,copyThreshHigh,copyThreshLow,PDmode,namesFreq,binSize,lambda,age,lifetimes,age_count)

nBins = nPop/binSize;

% death - birth
nBirth = binornd(nPop,pDeath);
indexDeath = randsample(nPop,nBirth);

% distribute births across bins
nBirthBinV = repmat(floor(nBirth/nBins),1,nBins);
rest = mod(nBirth,nBins);
if rest > 0
    binIndices = randperm(nBins,rest);
    nBirthBinV(binIndices) = nBirthBinV(binIndices)+1;
end

% divide population into bins
binM = reshape(randperm(nPop),binSize,nBins)';

% define copy pool
if copyAll == 0
    copyMask = (pop(2,:)>(t-copyThreshHigh)) & (pop(2,:)<(t-copyThreshLow));
    i = 1;
    while ~any(copyMask)
        copyMask = (pop(2,:)>(t-(copyThreshHigh+i))) & (pop(2,:)<(t-(copyThreshLow-i)));
        i = i+1;
    end
    copyIndex = find(copyMask);
else
    copyIndex = 1:nPop;
end

types = unique(pop(1,copyIndex)); % unique variant types in the copy pool
if numel(types)>1
    [~, ~, ic] = unique(pop(1,copyIndex));
    h = accumarray(ic,1)';
    if lambda == 0
        h = (h./numel(copyIndex));
    else
        [~, idx] = ismember(types,age(1,1:age_count));
        h_age = t-age(2,idx);
        h = (h./numel(copyIndex)).*(1-exp(-lambda*h_age));
        h = h./sum(h);
    end

    % Pre-allocate space
    hAdd = zeros(1,nBirth);
    hAddIdx = 0;
    hAddPlus = 0;

    % Sample for each bin
    for i = 1:nBins
        mask = unique(pop(1,binM(i,:)));
        % Sample variants
        %samples = randsrc(nBirthBinV(i),1,[types;h])';
        samples = randsample(types,nBirthBinV(i),true,h);

        % remove forbidden types
        forbidden = ismember(samples,mask);
        validSamples = samples(~forbidden);

        % collect valid copies
        nValid = numel(validSamples);
        hAdd(hAddIdx+1:hAddIdx+nValid) = validSamples;
        hAddIdx = hAddIdx + nValid;

        % count innovations needed
        hAddPlus = hAddPlus+sum(forbidden);
    end
    hAdd = hAdd(1:hAddIdx);
else
    hAddPlus = nBirth;
    hAdd = [];
end

% update population
newVariants = value + (1:hAddPlus);
pop(1,indexDeath) = [hAdd,newVariants];
pop(2,indexDeath) = t;

% update age tracking
age(1, age_count+1:age_count+hAddPlus) = newVariants;
age(2, age_count+1:age_count+hAddPlus) = t;
age_count = age_count+hAddPlus;

% check for extinctions and record lifetime
types = unique(pop(1,:));
[~, idx] = ismember(age(1,1:age_count),types);
dead_idx = find(idx==0);

if ~isempty(dead_idx)
    dead_variants = age(:,dead_idx);
    h_lifetimes = t-dead_variants(2,:);

    % emergency expansion of lifetimes-vector
    max_lifetime = max(h_lifetimes);
    if max_lifetime>length(lifetimes)
        new_size = max(length(lifetimes)*2,max_lifetime);
        lifetimes = [lifetimes,zeros(1,new_size-length(lifetimes))];
%        warning('lifetimes expanded to %d at t=%d', new_size, t);
    end
    % update lifetimes
    lifetimes = lifetimes+accumarray(h_lifetimes',1,[length(lifetimes),1])';

    % remove dead variants from age array
    age(:,dead_idx) = [];
    age_count = age_count-numel(dead_idx);
end


% update progeny frequency
if PDmode == 1
    names = unique(hAdd);
    [~, ~, ic] = unique(hAdd);
    progFreq = accumarray(ic, 1);
    namesFreq(names) = namesFreq(names)+progFreq';

    % Add innovations to progeny count
    namesFreq(newVariants) = 1;
end

value = value+hAddPlus;
addMut = hAddPlus/nBirth;

