%% SET FILE NAMES

% Create excel file names

excel    = [caseName '.xlsx']; 

data_sg  = [caseName '_data_sg.xlsx'];
data_vsc = [caseName '_data_vsc.xlsx'];
data_ipc = [caseName '_data_ipc.xlsx'];
data_shunt = [caseName '_data_shunt.xlsx'];

[folder_org, name_org, ext] = fileparts(which(excel));

data_sg2 = fullfile(folder_org, [name_org '_data_sg' ext]);
data_vsc2 = fullfile(folder_org, [name_org '_data_vsc' ext]);
data_ipc2 = fullfile(folder_org, [name_org '_data_ipc' ext]);
data_shunt2 = fullfile(folder_org, [name_org '_data_shunt' ext]);



% Simulink models names
linear    =  [caseName '_LIN'];
nonlinear = [caseName];

% Flag to indicate if T_case should be used (Obsolete)
shared_power = 0;