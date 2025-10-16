function [x, cost_reg_history, cost_data_history] =  delta2chi_MEDI(deltaB, Params, m, Mask, lambda, merit, edgePer)  
% [x, cost_reg_history, cost_data_history] =  delta2chi_MEDI(deltaB, Params, m, maskErode, lambda, merit, edgePer)   
%   QSM using MEDI_L1
%   deltaB: field change in ppm
%   Params: parameters structures
%   m: data term weighting
%   maskErode: brain Mask
%   lambda: for regularization
%   merit: fine tuning
%   edgePer: percentage of voxel as structure edges
% 
% modified according to MEDI_L1.m by Tian Liu, Ref: T. Liu et al. MRM 2013;69(2):467-76
% http://weill.cornell.edu/mri/pages/qsm.html
%%%%%%%%%%%%%%% weights definition %%%%%%%%%%%%%%
if nargin < 7
    edgePer = 0.9;
elseif nargin < 6
    merit = 0;
    edgePer = 0.9;    
elseif nargin < 5
    lambda = 1000;
    merit = 0;
    edgePer = 0.9;
end

% internal parameters for L1 solver
cg_max_iter = 100;              % for cg loop
cg_tol = 0.01;                  % cg tolerance

max_iter = 10;                  % outer loop
tol_norm_ratio = 0.1;

gradient_weighting_mode = 1;

grad = @cgrad;
div = @cdiv;

% kernel
D = conv_kernel_rot_c0(Params, Params.TAng);     
D = ifftshift(D);

% b0 and weighting
b0 = m.*exp(1i*deltaB);                         % change RDF to deltaB in ppm
wG = gradient_mask(gradient_weighting_mode, m, Mask, grad, Params.voxSize, edgePer);

iter=0;
x = zeros(Params.sizeVol); 
res_norm_ratio = Inf;
cost_data_history = zeros(1,max_iter);
cost_reg_history = zeros(1,max_iter);

e=0.000001; %a very small number to avoid /0

% main loop
while (res_norm_ratio>tol_norm_ratio)&&(iter<max_iter)
tic
    iter=iter+1;
    Vr = 1./sqrt(abs(wG.*grad(real(x),Params.voxSize)).^2+e);       % 1/abs(M*G*x_n)
    
    % complex
    w = m.*exp(1i*real(ifftn(D.*fftn(x))));                           % w: W*exp(i*D*x_n)
    
    % regularization term
    reg = @(dx) div(wG.*(Vr.*(wG.*grad(real(dx),Params.voxSize))),Params.voxSize);  
    
    % data fidelity term
    fidelity = @(dx)2*lambda*real(ifftn(D.*fftn(conj(w).*w.*real(ifftn(D.*fftn(dx))))));

    A =  @(dx) reg(dx) + fidelity(dx);       
    b = reg(x) + 2*lambda*real(ifftn(D.*fftn( conj(w).*conj(1i).*(w-b0))));

    dx = real(cgsolve(A, -b, cg_tol, cg_max_iter, 10));      % CG solver, solve A*dx = b
    res_norm_ratio = norm(dx(:))/norm(x(:));
    x = x + dx;

    wres=m.*exp(1i*(real(ifftn(D.*fftn(x))))) - b0;         % weighted residual, should be small

    cost_data_history(iter) = norm(wres(:),2);
    cost=abs(wG.*grad(x));
    cost_reg_history(iter) = sum(cost(:));

    
    if merit
        wres = wres - mean(wres(Mask(:)==1));   % unbias
        a = wres(Mask(:)==1);
        factor = std(abs(a))*6;     % 6 sigma
        
        wres = abs(wres)/factor;    % normalization        
        wres(wres<1) = 1;
        
        m = m./(wres.^2);        
        b0 = m.*exp(1i*deltaB);
        
    end
    
    fprintf('iter: %d; res_norm_ratio:%8.4f; cost_L2:%8.4f; cost:%8.4f.\n',iter, res_norm_ratio,cost_data_history(iter), cost_reg_history(iter));
toc
    
end

x = x.*Mask;

end

