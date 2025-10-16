% Discrete Divergence Using Backward Difference
% with the Dirichlet Boundary Condition

% Created by Youngwook Kee (Oct 21 2015)
% Last modified date: Oct 24 2015

% References:
% [1] Chambolle, An Algorithm for Total Variation Minimization and
% Applications, JMIV 2004

function div = bdiv3(Gx, voxel_size)

if (nargin < 2)
    voxel_size = [1 1 1];
end

sz=size(Gx);
ndim=numel(sz);
crd0=repmat({':'},[1 ndim]);
pd0=zeros([1 ndim]);
div=zeros([sz(1:3) 1 sz(5:end)]);
for i=1:3
    crd=crd0;crd{i}=1:sz(i)-1;crd{4}=i;
    pd=pd0;pd(i)=1;
    G=Gx(crd{:});
    div = div - (padarray(G,pd,0,'post')-padarray(G,pd,0,'pre'))/voxel_size(i);
end
div=reshape(div,[sz(1:3) sz(5:end)]);

end
