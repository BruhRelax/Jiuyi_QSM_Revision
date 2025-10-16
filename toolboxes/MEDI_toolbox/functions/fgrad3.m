% Discrete Gradient Using Forward Differences
% with the Neuman Boundary Condition

% Created by Youngwook Kee (Oct 21 2015)
% Last modified date: Oct 24 2015

% References:
% [1] Chambolle, An Algorithm for Total Variation Minimization and
% Applications, JMIV 2004
% [2] Pock et al., Global Solutions of Variational Models with Convex
% Regularization, SIIMS 2010

function Gx = fgrad3(chi, voxel_size)

if (nargin < 2)
    voxel_size = [1 1 1];
end

% chi = double(chi);
sz=size(chi);
ndim=numel(sz);
crd0=repmat({':'},[1 ndim]);
szA=[sz(1:3) 1 sz(4:end)];
for i=1:3
    crd=crd0;crd{i}=[2:sz(i) sz(i)];
    A{i}=reshape((chi(crd{:})-chi)/voxel_size(i),szA);
end
Gx=cat(4,A{:});
end
