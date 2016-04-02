function [nevList] = findMSIT(patientID,patientDirectory)
%FINDMSIT finds MSIT data in a patient directory
%   nevList = findMSIT(patientID,patientDirectory) will search recursively through the directory 
% 	in patientDirectory for MSIT data and will return the names of files that have
%	the trigger structure output from psychToolbox. These nev files will also be copied to 
% 	a directory in the shethLabBackup Data folder called [patientID]_MSIT.

% author: EHS20160210  

cd /mnt/mfs/shethLabBackup_nov15/Code/Analysis/

% 1) find .nev files in the directory
dirlist = subdir(fullfile(patientDirectory,'*.nev'));

% 2) open NEV files
addpath(genpath('./NPMK'))

% initializing session count and list of files
sessionCount = 0;
nevList = cell(0);
nevListStr = char();

% looping over nev files.
for fl = 1:length(dirlist)
    NEV = openNEV(dirlist(fl).name,'nomat','nosave');
    triggers = NEV.Data.SerialDigitalIO.UnparsedData;
    if isempty(triggers)
        display(sprintf('no triggers found in file %s',dirlist(fl).name))
    elseif isequal(triggers(1),255) && isequal(triggers(2),90)
        sessionCount = sessionCount+1;
        display('found the start of an MSIT session!')
        N = sum(triggers==90);
        display(sprintf('found %d trials of MSIT in file %s',N,dirlist(fl).name))
        nevList = vertcat(nevList,dirlist(fl).name);
	nevListStr = horzcat(nevListStr,' ',dirlist(fl).name);
    elseif isequal(triggers(1),255) && ~isequal(triggers(2),90)
	display('There is the start of a task in this file, though it isn't MSIT.')
    else
        display(sprintf('There are triggers in file %s, though they might not be MSIT.',dirlist(fl).name))
    end
end


eval(['! mkdir /mnt/mfs/shethLabBackup_nov15/Data/' patientID '/'])
eval(['! cp -v ' nevListStr '  /mnt/mfs/shethLabBackup_nov15/Data/' patientID '/'])

end
