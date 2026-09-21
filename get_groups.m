function [sizes, bnd, perm, g] = get_groups(nPop, groupSize)

x = max(1, round(nPop/groupSize));   % closest group count to target
q = floor(nPop/x);                   % base size
r = mod(nPop,x);                     % this many groups get one extra

sizes = [repmat(q+1,1,r), repmat(q,1,x-r)];
bnd   = [0 cumsum(sizes)];

% assign individuals to groups
perm = randperm(nPop);

if nargout > 3
    g = zeros(1,nPop);
    g(perm) = repelem(1:numel(sizes), sizes);
end
