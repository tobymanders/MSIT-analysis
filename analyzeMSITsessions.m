function [sessionInfo] = analyzeMSITsessions(patientID)
% ANALYZEMSITSESSIONS doesg MSIT analysis over sessions.
%   Finds and enters directory called [PatientID].
%   Looks for data in subdirectory ./Data and analyzes that data.
%   MSIT over behavioral sessions.


if ~exist(patientID,'dir')
    display('no patient directory found')
else
    % entering patient directory
    fullPath = pwd;
    dirLen = length(patientID);
    if ~strcmp(fullPath(end-(dirLen-1):end),patientID)
        cd(patientID)
        fullPath = pwd;
    end
    
    % creating file structure
    if ~exist('Data','dir')
        display(sprintf('making data directory. You should probably put your data in %s',strcat(fullPath,'Data')))
        mkdir(fullPath,'Data');
        if isunix
            addpath('./Data')
        else
            addpath('.\Data')
        end
    end
    if ~exist('Figs','dir')
        mkdir(fullPath,'Figs');
    end
    
    % finding .nev files
    dirList = dir(fullfile(fullPath,'Data'));
    display(['found ' num2str(length(dirList)-2) ' MSIT sessions!'])
    
    % looping over sessions
    for  ss = 1:length(dirList)-2
        nevFile = fullfile(fullPath,'Data',dirList(ss+2).name);
        
        % analyzing behavior and saving figures and stats
        try
            [conflictStats{ss},FBstats{ss}]...
                = analyzeMSITbehavior(patientID,ss,nevFile);
        catch
            display(sprintf('problem analyzing data from session %s.', ss))
        end
        
        % analyzing neurons and saving figures and stats
        try
            [neuralStats{ss}] = analyzeMSITunits (patientID,ss,nevFile);
        catch
            display(sprintf('problem analyzing data from session %s.', ss))
        end
        
        % saving statistics.
        saveFlag = 1;
        if saveFlag
            save(['./Data/' patientID '_statisticsOverSessions.mat'],'conflictStats','FBstats','neuralStats')
        end
        
    end
end

% create summary report

end

