% ANALYZEMSIT will analyze all of the MSIT data in a folder, provided that
% the data is organized correctly.
%
% How to organize:
%   1) make a folder for each subject named by their identifier.
%   2) put the data from each subject in a folder called 'Data' in each of
%       corresponding patient's folders.
%
% Running this or any of the other 'MSIT-analysis' scripts will also create
% a folder called 'Figs' in each subject's folder which will contain
% figures of analyses for that specific subject.
%
% statistics will be saved in each subject's 'Data' folder and will be
% visualized in their figures.
%
%


% author: ElliotHSmith
% https://github.com/elliothsmith/MSIT-analysis


folder_name = uigetdir('/home/patients','Which directory contains subject folders?');

ptlist = dir(folder_name);
for pt = 1:length(ptlist.name)
    
    [ptPath, ptName, nvExt] = fileparts(ptlist(pt).name);
    
    sessionList = dir(['/home/elliot/data/msit_units/' ptName '/Data/*.nev']);
    for ss = 1:length(sessionList.name)
        
%         [behavioralStats] = analyzeMSITbehavior(patientID,session,nevFile)
%         
%         [neuralStats] = analyzeMSITunits (patientID,session,nevFile)
        
        coherenceStats{pt,ss} = analyzeMSITcoherence(ptName,ss,nevFile)
        
    end
    
    
    
    
end


