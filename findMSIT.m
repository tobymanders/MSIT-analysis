function [nevList] = findMSIT(patientID,patientDirectory,nsFlag)
%FINDMSIT finds MSIT data in a patient directory
%   nevList = findMSIT(patientID,patientDirectory) will search recursively through the directory
% 	in patientDirectory for MSIT data and will return the names of files that have
%	the trigger structure output from psychToolbox. These nev files will also be copied to
% 	a directory in the shethLabBackup Data folder called [patientID]_MSIT.
%
%   This code will also transfer the associated NS3 file for local field
%   potential analysis if teh nsFlag is set to 1 (default, 0).

% versionDate: 20160402
% author: EHS

cd /mnt/mfs/shethLab/Code/Analysis/

% setting defaults.
if nargin <3
    nsFlag = false;
end
logical nsFlag

% 1) find .nev files in the directory
dirlist = subdir(fullfile(patientDirectory,'*.nev'));

% 2) open NEV files
addpath(genpath('./NPMK'))

% initializing session count and list of files
sessionCount = 0;
nevList = cell(0);
nevListStr = char();
if nsFlag
    nsList = cell(0);
    nsListStr = char();
end

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
        if nsFlag
            nsList = vertcat(nsList,[dirlist(fl).name(1:end-3) '.ns3']);
            nsListStr = horzcat(nevListStr,' ',[dirlist(fl).name(1:end-3) '.ns3']);
        end
    elseif isequal(triggers(1),255) && ~isequal(triggers(2),90)
        display('There is the start of a task in this file, though it isn't MSIT.')
    else
        display(sprintf('There are triggers in file %s, though they might not be MSIT.',dirlist(fl).name))
    end
end

display(sprintf('copying NEV files to ~/Data/%s/',patientID))
eval(['! mkdir ~/Data/' patientID '/'])
eval(['! cp -v ' nevListStr '  ~/Data/' patientID '/'])

if nsFlag
    display(sprintf('copying NS3 files to ~/Data/%s/',patientID))
    eval(['! mkdir ~/Data/' patientID '/'])
    eval(['! cp -v ' nsListStr '  ~/Data/' patientID '/'])
end

return


function [sub,fls] = subdir(CurrPath)
%   SUBDIR  lists (recursive) all subfolders and files under given folder
%
%   SUBDIR
%        returns all subfolder under current path.
%
%   P = SUBDIR('directory_name')
%       stores all subfolders under given directory into a variable 'P'
%
%   [P F] = SUBDIR('directory_name')
%       stores all subfolders under given directory into a
%       variable 'P' and all filenames into a variable 'F'.
%       use sort([F{:}]) to get sorted list of all filenames.
%
%   See also DIR, CD

%   author:  Elmar Tarajan [Elmar.Tarajan@Mathworks.de]
%   version: 2.0
%   date:    07-Dez-2004
%
if nargin == 0
    CurrPath = cd;
end% if

if nargout == 1
    sub = subfolder(CurrPath,'');
else
    [sub fls] = subfolder(CurrPath,'','');
end% if


function [sub,fls] = subfolder(CurrPath,sub,fls)
%------------------------------------------------
tmp = dir(CurrPath);
tmp = tmp(~ismember({tmp.name},{'.' '..'}));
for i = {tmp([tmp.isdir]).name}
    sub{end+1} = [CurrPath '\' i{:}];
    if nargin==2
        sub = subfolder(sub{end},sub);
    else
        tmp = dir(sub{end});
        fls{end+1} = {tmp(~[tmp.isdir]).name};
        [sub fls] = subfolder(sub{end},sub,fls);
    end% if
end% if
