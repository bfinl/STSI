Subject = subjID;

if(ispc)
% Current code is configured based on PC environment, 
% Please change the MainPath to where the folder 'STSI_Codes' is located.
    MainPath = 'E:\STSI_CodesAndNotes\STSI_Codes\';

    % These toolboxes needs to be in the \ToolBox_External folder

    % FieldTrip
    cd(fullfile(MainPath, 'ToolBox_External'));
    folder_name_fieldtrip = dir(fullfile(pwd, '*fieldtrip*'));
    if isempty(folder_name_fieldtrip)
        error('FieldTrip folder not found in ToolBox_External. Please download and place it there.');
    end
    path_fieldtrip = fullfile(folder_name_fieldtrip(1).folder, folder_name_fieldtrip(1).name);
    addpath(path_fieldtrip);
    ft_defaults;

    % EEGLAB
    cd(fullfile(MainPath, 'ToolBox_External'));
    folder_name_eeglab = dir(fullfile(pwd, '*eeglab*'));
    if isempty(folder_name_eeglab)
        error('EEGLAB folder not found in ToolBox_External. Please download and place it there.');
    end
    path_eeglab = fullfile(folder_name_eeglab(1).folder, folder_name_eeglab(1).name);
    cd(path_eeglab);
    eeglab;

    % TensorLab
    cd(fullfile(MainPath, 'ToolBox_External'));
    folder_name_tensorlab = dir(fullfile(pwd, '*tensorlab*'));
    if isempty(folder_name_tensorlab)
        error('Tensorlab folder not found in ToolBox_External. Please download and place it there.');
    end
    path_tensorlab = fullfile(folder_name_tensorlab(1).folder, folder_name_tensorlab(1).name);
    addpath(path_tensorlab)


    Code_Folder = [MainPath '\Codes'];
    addpath(genpath(Code_Folder))

    PreProcess_Folder = [MainPath '\Preprocessing\'];

    GridLoc_Folder = [MainPath '\Grid_Location_Parameters\' Subject '\'];
    if ~exist(GridLoc_Folder, 'dir')
       mkdir(GridLoc_Folder)
    end
    
    Fig_Folder = [MainPath '\Figures\Spike_',Subject '\'];
    if ~exist(Fig_Folder, 'dir')
       mkdir(Fig_Folder)
    end

    Spk_Folder = [MainPath '\Spikes\Spike_',Subject '\'];
    if ~exist(Spk_Folder, 'dir')
       mkdir(Spk_Folder)
    end

    Sz_Folder = [MainPath '\Seizures\Seizure_',Subject '\'];
    if ~exist(Sz_Folder, 'dir')
      mkdir(Sz_Folder)
    end

    savepath;
end
