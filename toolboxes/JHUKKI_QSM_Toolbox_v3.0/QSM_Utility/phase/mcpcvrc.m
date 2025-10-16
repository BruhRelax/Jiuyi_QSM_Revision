function GREPhase = mcpcvrc(GREMag, GREPhase, Params)
% GREPhase = mcpcvrc(GREMag, GREPhase, Params)
% 
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% Multiple channel phase combination - using virtual reference coil method
% (vrc method), using mcpcc to get the virtual reference first
% Ref. Robinson et al, NMRB, 2016
%        Parker et al, MRM, 2014
% updated 2017-03-15, X.L.

%% get virtual reference phase using mcpcc method first
centVox = floor(0.5*(Params.sizeVol));    
Nc = Params.nchannel;
Ne = Params.nEchoes;

if Params.nEchoes > 1   % multiecho data X*Y*Z*NE*Nc    
    % do global phaes correction, reference to center voxel    
    for IndEcho = 1:Ne
        for IndChannel = 1:Nc

            GREPhase(:,:,:,IndEcho,IndChannel) = ...
            angle(exp(1i.*GREPhase(:,:,:,IndEcho,IndChannel)).*...
                (exp(-1i.*GREPhase(centVox(1), centVox(2), centVox(3), IndEcho, IndChannel))));    
                
        end
    end
    
    VRC = sum(GREMag.*exp(1i.*GREPhase), 5)./sum(GREMag, 5);    % virtual reference coil
    VRC(isnan(VRC)) = 0;
           
    PhaseDiff = angle(exp(1i.*GREPhase).*repmat(conj(VRC), [1, 1, 1, 1, Nc])); % coil phase difference
    
else                    % single echo data
    % reference to a voxel with good SNR for every coil
    GREMagCoilMin = min(GREMag, [], 4);
    [~, Ind] = max(GREMagCoilMin(:));
    [centVox(1), centVox(2), centVox(3)] = ind2sub(Params.sizeVol, Ind);
    
    % do global phaes correction, reference to center voxel    
    for IndChannel = 1:Nc
        GREPhase(:,:,:,IndChannel) = ...
        angle(exp(1i.*GREPhase(:,:,:,IndChannel)).*(exp(-1i.*GREPhase(centVox(1), centVox(2), centVox(3), IndChannel))));
    end
    
    VRC = sum(GREMag.*exp(1i.*GREPhase), 4)./sum(GREMag, 4);    % virtual reference coil
    VRC(isnan(VRC)) = 0;
    
    PhaseDiff = angle(exp(1i.*GREPhase).*repmat(conj(VRC), [1, 1, 1, Nc])); % coil sensitivity phase difference
end

%% get phase difference between each coil and VRC and do lowpass filter 
if Params.nEchoes > 1
    for IndEcho = 1:Ne
       for IndChannel = 1:Nc
            PhaseDiff(:,:,:,IndEcho,IndChannel) = smooth3(PhaseDiff(:,:,:,IndEcho,IndChannel));     % smooth3 or medfilt3
       end
    end
    GREPhase = angle(exp(1i.*GREPhase).*exp(-1i.*PhaseDiff));        % normalized to VRC    
    GREPhase = angle(sum(GREMag.*exp(1i.*GREPhase), 5));        % multiple coils combined
        
else
    
    for IndChannel = 1:Nc
        PhaseDiff(:,:,:,IndChannel) = smooth3(PhaseDiff(:,:,:,IndChannel));   
    end
    GREPhase = angle(exp(1i.*GREPhase).*exp(-1i.*PhaseDiff));         % normalized to VRC    
    GREPhase = angle(sum(GREMag.*exp(1i.*GREPhase), 4));              % multiple coils combined 
end



