function UpdateTable(handles, statustext )
    %% Author: Jiri van Bergen
    % Affiliation: Radiology @ JHU - Kirby Center
    % Contact via xuli@mri.jhu.edu
    
    % Update table with new status for this dataset
    curTableData = get(handles.TableDatasets, 'Data');
    curTableData{handles.CurrentDataset, 2} = statustext;
    set(handles.TableDatasets, 'Data', curTableData);
end

