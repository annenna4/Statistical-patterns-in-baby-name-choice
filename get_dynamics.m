function [pop,value,namesFreq,age,lifetimes,age_count] = get_dynamics_novelty2(t,pop,value,pDeath,nPop,pMut,lambda,copyAll,copyThreshHigh,copyThreshLow,PDmode,namesFreq,age,lifetimes,age_count)

% death
nBirth = binornd(nPop,pDeath); % number of death = number of birth
if nBirth == 0
    return; % early exit if no births
end
indexDeath = randsample(nPop,nBirth); % generating indices of individuals to be removed

% reproduction
nMut = binornd(nBirth,pMut); % number of innovations
nCopy = nBirth-nMut; % number of copies

% emergency expansion of age-vector and namesFreq-vector
needed_age_cols = age_count+nMut;
if needed_age_cols > size(age, 2)
    new_size = max(size(age,2)*2,needed_age_cols);
    age = [age,zeros(2,new_size-size(age,2))];
    warning('age expanded to %d columns at t=%d',new_size,t);
end
if PDmode == 1
    needed_freq_size = value + nMut;
    if needed_freq_size > length(namesFreq)
        new_size = max(length(namesFreq)*2,needed_freq_size);
        namesFreq = [namesFreq,zeros(1,new_size-length(namesFreq))];
        warning('namesFreq expanded to %d at t=%d',new_size,t);
    end
end

if nCopy > 0
    if copyAll == 0
        copyMask = (pop(2,:)>(t-copyThreshHigh)) & (pop(2,:)<(t-copyThreshLow));
        i = 1;
        while ~any(copyMask)
            copyMask = (pop(2,:)>(t-(copyThreshHigh+i))) & (pop(2,:)<(t-(copyThreshLow-i)));
            i = i + 1;
        end
        copyIndex = find(copyMask);
    else
        copyIndex = 1:nPop;
    end

    if lambda==0 % unbiased transmission
        h = numel(copyIndex);
        index = randsample(h,nCopy,true);
        hAdd = pop(1,copyIndex(index));
    else % anti-novelty bias of strength lambda
        % choosing role models from copy pool
        types = unique(pop(1,copyIndex)); % variant types present in copy pool
        if numel(types)>1
            [~, ~, ic] = unique(pop(1,copyIndex));
            h = accumarray(ic,1)';

            [~, idx] = ismember(types,age(1,1:age_count));
            h_age = t-age(2,idx);
            h = (h./numel(copyIndex)).*(1-exp(-lambda*h_age));
            hAdd = randsrc(1,nCopy,[types;h./sum(h)]);
        else
            hAdd = types*ones(1,nCopy);
        end
    end
else
    hAdd = [];
end

% update populations
if nMut > 0
    new_variants = value+(1:nMut);
    pop(1,indexDeath) = [hAdd,new_variants];
else
    pop(1,indexDeath) = hAdd;
end
pop(2,indexDeath) = t; % adding birth years

% update age
if nMut > 0
    age(1, age_count+1:age_count+nMut) = new_variants;
    age(2, age_count+1:age_count+nMut) = t;
    age_count = age_count + nMut;
end

% check for extinctions and record lifetime
types = unique(pop(1,:));
[~, idx] = ismember(age(1,1:age_count),types);
dead_idx = find(idx == 0);

if ~isempty(dead_idx)
    dead_variants = age(:,dead_idx);
    h_lifetimes = t - dead_variants(2,:);

    % emergency expansion of lifetimes-vector
    max_lifetime = max(h_lifetimes);
    if max_lifetime > length(lifetimes)
        new_size = max(length(lifetimes)*2,max_lifetime);
        lifetimes = [lifetimes, zeros(1,new_size-length(lifetimes))];
        warning('lifetimes expanded to %d at t=%d',new_size,t);
    end
    % update lifetimes
    lifetimes = lifetimes+accumarray(h_lifetimes',1,[length(lifetimes),1])';

    % remove dead variants from age array
    age(:, dead_idx) = [];
    age_count = age_count-numel(dead_idx);
end

% update progeny frequency
if PDmode == 1 && ~isempty(hAdd)
    names = unique(hAdd);
    [~, ~, ic] = unique(hAdd);
    progFreq = accumarray(ic,1);
    namesFreq(names) = namesFreq(names)+progFreq';

    % Add innovations to progeny count
    if nMut > 0
        namesFreq(value+(1:nMut)) = 1;
    end
end

value = value+nMut; % update variant count
