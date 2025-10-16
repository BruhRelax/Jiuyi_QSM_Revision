function [unwrappedPhase]  = phase_unwrap_laplacian(dataPhase, Params, RefVox, zp, mask)
% [unwrappedPhase]  = phase_unwrap_laplacian(dataPhase, Params, RefVox, zp, mask)
%================================
% HEADER
%================================
%
% Name:             phase_unwrap_laplacian.m
% Author:           Xu Li, PhD
% Email:            xuli@mri.jhu.edu
%
%================================
% PURPOSE
% unwrap phase according to the laplacian method
% i.e. del2(theta) = cos(theta)*del2(sin(theta)) + sin(theta)*del2(cos(theta))
%  Updated by Xu Li, 2019-06-20

textWaitbar = 'Performing phase unwrapping';
multiWaitbar(textWaitbar, 0, 'Color', 'b' );

warning off all

%% zero padding if necessary
N = size(dataPhase); 
Ny = N(1);
Nx = N(2);
Nz = N(3);

if nargin < 3      
    zp = 1;         % 1 for new code and 2 for old 
    RefVox = [floor(N(1)/2), floor(N(2)/2), floor(N(3)/2)]; 
    mask = ones(size(dataPhase));
elseif nargin < 4
    zp = 1;    
    mask = ones(size(dataPhase));
elseif nargin < 5
    mask = ones(size(dataPhase));    
end

if isempty(RefVox)
    RefVox = [floor(N(1)/2), floor(N(2)/2), floor(N(3)/2)]; 
end

RefFlag = 1;                % Do Ref in general
if RefVox == 0         
    RefFlag = 0;            % do not do Ref
end

if zp == 1
    % check if there are odd dimensions and do padding
    % only slice number could be odd number
    
    if max(N) < 128
        padsize = [32,32,32];
    elseif max(N) < 256     % higher resolution
        padsize = [64,64,64];        
    else % ultra-high        
        padsize = [128, 128, 128];
    end
    
    dimsOdd = mod(size(dataPhase), 2);
    
    dataPhase = padarray(dataPhase, dimsOdd, 'replicate', 'post');  % padding to even dim first     
    dataPhase = padarray(dataPhase, padsize, 'replicate');
    
    mask = padarray(mask, dimsOdd, 'replicate', 'post');  % padding to even dim first     
    mask = padarray(mask, padsize, 'replicate');    
    
    N = size(dataPhase);                        % new dimentions after padding, now all even

    ksize = [3, 3, 3];               
    khsize = (ksize-1)/2;
    
    kernelcore = [];
    kernelcore(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
    kernelcore(:,:,2) = [0 1 0; 1 -6 1; 0 1 0];
    kernelcore(:,:,3) = [0 0 0; 0 1 0; 0 0 0];

    Kernel = zeros(N);
    Kernel( 1+N(1)/2 - khsize(1) : 1+N(1)/2 + khsize(1), 1+N(2)/2 - khsize(2) : 1+N(2)/2 + khsize(2), ...
        1+N(3)/2 - khsize(3) : 1+N(3)/2 + khsize(3) ) = kernelcore;

    % laplacian kernel in k space
    lap_opK = fftn(fftshift(Kernel));                                     % make real delta function, only do fftn

    multiWaitbar(textWaitbar, 0.2, 'Color', 'b' );
    
    % inverse laplacian kernel       
    inv_lap_opK = zeros(size(lap_opK));
    inv_lap_opK(lap_opK~=0) = 1./lap_opK(lap_opK~=0);                   % inverse of kernel in k-space

    multiWaitbar(textWaitbar, 0.4, 'Color', 'b' );
    
    % laplacian based phase unwrapping 
    PhaseSin = sin(dataPhase);
    PhaseCos = cos(dataPhase);    
    PhaseLaplacian  = PhaseCos.*ifftn(  lap_opK.* fftn(PhaseSin)) - PhaseSin.*ifftn( lap_opK.* fftn(PhaseCos));  %  
    
    clear PhaseSin PhaseCos
    multiWaitbar(textWaitbar, 0.8, 'Color', 'b' );
    
    unwrappedPhase  =  ifftn( inv_lap_opK.* fftn(PhaseLaplacian.*mask));   %      
    unwrappedPhase = unwrappedPhase(padsize(1)+(1:Ny), padsize(2)+(1:Nx), padsize(3)+(1:Nz));               % crop back
    
else
    % old code, obsolete
    % and calculate laplacian in k-space using k-space cooridnate
    Nx = (Params.sizeVol(2));
    Ny = (Params.sizeVol(1));
    Nz = (Params.sizeVol(3));

    FOVx = Params.fov(2);
    FOVy = Params.fov(1);
    FOVz = Params.fov(3);

    Nx_zp = Nx*zp;
    Ny_zp = Ny*zp;
    Nz_zp = Nz*zp;
    
    zpx = Nx_zp/Nx;                     % final zero padding ratio
    zpy = Ny_zp/Ny;
    zpz = Nz_zp/Nz;

    FOVx_zp = FOVx*zpx;
    FOVy_zp = FOVy*zpy;
    FOVz_zp = FOVz*zpz;

    x1 = floor((Nx_zp - Nx)/2) + 1;     % index where to put the original data
    x2 = x1 + Nx - 1;
    y1 = floor((Ny_zp - Ny)/2) + 1;
    y2 = y1 + Ny - 1;
    z1 = floor((Nz_zp - Nz)/2) + 1;
    z2 = z1 + Nz - 1;

    Phasej = dataPhase;
    if isa(dataPhase(1), 'single')
        dataPhase = zeros(Ny_zp, Nx_zp, Nz_zp, 'single');
    else
        dataPhase = zeros(Ny_zp, Nx_zp, Nz_zp);    
    end
    % put the data in the center of the zero padded one
    dataPhase(y1:y2, x1:x2, z1:z2) = Phasej;
    clear Phasej
    
    %-------------------   calculate the k^2 in k space
    % Get k space cooridnate
    dkx = 1/FOVx_zp;
    dky = 1/FOVy_zp;
    dkz = 1/FOVz_zp;

    kx = linspace(-Nx_zp/2+1, Nx_zp/2, Nx_zp).*dkx;
    ky = linspace(-Ny_zp/2+1, Ny_zp/2, Ny_zp).*dky;
    kz = linspace(-Nz_zp/2+1, Nz_zp/2, Nz_zp).*dkz;
    
    % Update
    multiWaitbar(textWaitbar, 0.2, 'Color', 'b' );

    if isa(dataPhase(1), 'single')
        % Form the KSq_Grid using a for loop
        KSq_Grid = zeros(Ny_zp, Nx_zp, Nz_zp, 'single');
        for kk = 1:Nz_zp
            for jj = 1:Nx_zp
                for ii = 1:Ny_zp
                    KSq_Grid(ii, jj, kk) = (kx(jj)^2 + ky(ii)^2 + kz(kk)^2);
                end
            end
        end
    else
        [KX_Grid, KY_Grid, KZ_Grid] = meshgrid(kx, ky, kz);  % mesh in k space
        KSq_Grid = (KX_Grid.^2 + KY_Grid.^2 + KZ_Grid.^2);

        clear KX_Grid KY_Grid KZ_Grid
    end

    % Update
    multiWaitbar(textWaitbar, 0.45, 'Color', 'b' );

    %% calculate 
    ThetaSin = sin(dataPhase);
    ThetaCos = cos(dataPhase);
    clear dataPhase

    laplacianSin = fftn(ThetaSin);
    laplacianSin = fftshift(laplacianSin).*KSq_Grid;
    laplacianSin = ifftshift(laplacianSin);
    laplacianSin = ifftn(laplacianSin);
    laplacianSin = real(laplacianSin);

    % Update
    multiWaitbar(textWaitbar, 0.60, 'Color', 'b' );

    laplacianCos = fftn(ThetaCos);
    laplacianCos = fftshift(laplacianCos).*KSq_Grid;
    laplacianCos = ifftshift(laplacianCos);
    laplacianCos = ifftn(laplacianCos);
    laplacianCos = real(laplacianCos);

    % Update
    multiWaitbar(textWaitbar, 0.75, 'Color', 'b' );

    A = laplacianSin.*ThetaCos - laplacianCos.*ThetaSin;
    clear ThetaSin ThetaCos laplacianSin laplacianCos

    % Update
    multiWaitbar(textWaitbar, 0.85, 'Color', 'b' );

    unwrappedPhase = fftshift(fftn(A))./KSq_Grid;
    unwrappedPhase(isinf(unwrappedPhase)) = 0;
    clear A KSq_Grid

    % Update
    multiWaitbar(textWaitbar, 0.95, 'Color', 'b' );

    unwrappedPhase = real(ifftn(ifftshift(unwrappedPhase)));
    unwrappedPhase(isnan(unwrappedPhase)) = 0;
    
    % crop back to original size
    unwrappedPhase = unwrappedPhase(y1:y2, x1:x2, z1:z2);

    
end

%% shift to a known constant (find the largest phase region and select the middle point for reference)
                       % use shift in most cases
if RefFlag == 1
    c = - unwrappedPhase(RefVox(1), RefVox(2), RefVox(3));
else
    c = 0;
end

% Update
multiWaitbar(textWaitbar, 0.99, 'Color', 'b' );

unwrappedPhase = unwrappedPhase + c;        % phi' as in Schofield & Zhu, Optics Letters, 2003