function y = ss_mask(x, m)
if mod(numel(x),numel(m))~=0
    error(['Incompatible size of x ' num2str(size(x)) ' and mask ' num2str(size(m))]);
end
y=reshape(x, numel(m), []);
y=y(m>0,:);
end