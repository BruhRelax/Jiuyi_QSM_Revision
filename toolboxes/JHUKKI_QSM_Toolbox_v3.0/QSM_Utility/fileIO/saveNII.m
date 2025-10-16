function []=saveNII(saveVar, fileName, Params, permuteflag, ext, prec)
    % Wrapper function for NIFTI save operation
    % updated for v2.5, X.L., 2014-09-08
    % updated for v2.8, X.L., 2016-12-08
    % updated to include .img/.hdr format, X.L., 2017-05-31
    % 2019-02-18, to include dynamics data
    % 2019-06-05, included .nii.gz and prec, X.L.
    
    if nargin < 4
        permuteflag = 1;        % default setting, compatible with old version
        ext = '.nii.gz';
        prec = 16;
    elseif nargin < 5
        ext = '.nii.gz';           % default using nii
        prec = 16;
    elseif nargin < 6
        prec = 16;
    end
    
    if permuteflag == 1
        saveVar = permute(saveVar, [2, 1, 3:length(size(saveVar))]);
        Params.voxSize = Params.voxSize([2, 1, 3]);
    end
        
    if strcmpi(ext, '.nii') || strcmpi(ext, '.nii.gz')
        nii = make_nii(saveVar, Params.voxSize, [], prec);       % single, float 32
        % save as nifti
        nii.hdr.hist.qform_code = 1;
        nii.hdr.hist.quatern_d = 1;        
        save_nii(nii, [fileName, ext]);
    else
        % save as analyze
        nii = make_nii(flip(saveVar, 2), Params.voxSize, [], prec);
        save_nii(nii, [fileName, '']);
    end
    
    clear nii
end

