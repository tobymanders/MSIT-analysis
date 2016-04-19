function [coherenceStats] = analyzeMSITcoherence(patientID,sessionNum,nevFile)
%ANALYZEMSITCOHERENCE does spike-field coherence analysis on MSIT data.
%
%   [coherenceStats] = analyzeMSITcoherence(patientID,sessionNum,nevFile)
%   calculates spike-field coherograms for MSIT data in nevFile and its
%   associated ns(3) file. analyzeMSITcoherence saves statistics, and
%   plots results.
%
%   Two tips for easy usage:
%       1) place the nev file in the same directory as the nsx files.
%       2) make sure that the nev file with sorted units has its original
%           name at the beginning of the file name.
%


% author: ElliotHSmith (https://github.com/elliothsmith/MSIT-analysis)


%% loading data from NEV file
display('loading action potential and local field potential data...')
[dataPath, nvName, nvExt] = fileparts(nevFile);

% defining nsFile.
nsFile = fullfile(dataPath,[nvName(1:19) '.ns3']);

% parsing files and loading data.
if strcmp(nvExt,'.nev')
    NEV = openNEV(nevFile,'read');
    if ~exist(nsFile,'file')
        [nsFile,nsPath,~] = uigetfile('*.ns*','Select the correct NSx file');
        NS3 = openNSx(fullfile(nsPath,nsFile));
    else
        NS3 = openNSx(nsFile);
    end
elseif strcmp(nvExt,'.mat')
    load(nevFile);
end


%% organizing important task parameters.
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);


%% reaction time calculation
rt = trigTimes(trigs>=100 & trigs<105) - trigTimes(trigs>=1 & trigs<28);


%% defining Chronux parameters.
movingWin = [1 0.010];
params.Fs = 2e3; % sampling frequency for LFP
params.pad = 2; %
params.fpass = [0 50]; % frequency range of interest
params.tapers = [5 9]; % emphasize smoothing for the spikes
params.trialave = 0; % average over trials {CHANGES BELOW}
params.err = [2 0.01]; % population error bars


%% look for the start of the task
if isequal(trigs(1),255)
    display('task started in this recording')
else
    display('no task start found...')
    display('recording may have started after the behavioral task.')
end


%% timing (seconds)
pre = 2;
post = 3;


%% creating neural timing variable
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];
ChanUnitTimestamp(ChanUnitTimestamp(:,2)==255,:) = [];

inclChans = unique(ChanUnitTimestamp(:,1));
nChans = length(inclChans);
unitChanLabels = [NEV.ElectrodesInfo.ElectrodeLabel]';
unitChanLabels = strtrim(unitChanLabels(inclChans,:));


%% aligning data on specific or all trial markers: query user.
% alignSpot = input('\nAlign on: \n1)fixation? \n2)cue? \n3)response? \n4)all of the above ?');
% if isequal(alignSpot,4)
%     aS = 1:3;
% else
%     aS = alignSpot;
% end


display('Aligning AP data on stimulus and response.');
% for loop to save multiple epochs
for aS = 1:2
    % which alignment spot
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            nTrials = sum(trigs>=1 & trigs<28);
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>=100 & trigs<=105);
            nTrials = 300;
    end
    
    
    % looping over Channels
    for ch = 1:nChans
        
        % looping over number of units in the AP data
        nUnits = length(unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));
        for un = 1:nUnits
            
            % getting unit times for the current channel and unit.
            unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un,3); % in seconds
            
            % loooping over trials
            for tt = 1:nTrials
                
                %% putting the data in a structure.
                
                %%~~~~~~~~~~~~~~~~~~~~~~~~FROM O)LD CODE
                % data(ch).channel(un).unit(tt).times = unitTimes(unitTimes>trialStart-pre & unitTimes<trialStart+post) - repmat(trialStart,length(unitTimes(unitTimes>trialStart-pre & unitTimes<trialStart+post)),1);
                %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                
                %% TODO: Do I need to make sure that all of the times are positive?
                data(aS).channel(ch).unit(un).trial(tt).times = unitTimes(unitTimes>(trialStarts(tt)-pre) & unitTimes<trialStarts(tt)+post)...
                    - repmat(trialStarts(tt)-pre,length(sum(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post)),1);
                
                
            end
        end
    end
    
    % save data structure?
    saveFlag = 0;
    if saveFlag
        display(['saving spike time structure for ' alignName '...'])
        save([patientID '_session' num2str(sessionNum) '_spikeTimeStruct_alignedon' alignName '.mat'],'data','pre','post','alignName')
    end
    
end


%%  alinging LFP data on the same time base as AP data.
tmp = double(NS3.Data);
LFPlabels = {NS3.ElectrodesInfo.Label};
LFPmat = zeros(((pre+post)*params.Fs)+1,nTrials,length(LFPlabels),2);

%% TODO: determine LFP cahnnels based on labels and analyze proximal and distal trodes.
% BFtrodes = unique(unitChanLabels)
% proximalLFPtrode =  
% distalLFPtrode = 

%% TODO: figure out a smart way of allocating which LFP contacts you want
display('Aligning LFP data on stimulus and response.');
% for loop to save multiple epochs
for aS = 1:2
    % which alignment spot
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            nTrials = sum(trigs>=1 & trigs<28);
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>=100 & trigs<=105);
            nTrials = 300;
    end
    
    for ch = 1:length(LFPlabels)
        for tt = 1:nTrials
            
            %%~~~~~~~~~~~~~~~~~~~FROM OLD (WORKING) CODE
            % LFPmat(ch,:,trl) = ECoG2kp(ch,ceil(trialStart(trl)*Fs + win(1)*Fs):ceil(trialStart(trl)*Fs + win(2)*Fs));
            %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
            %% money step:
            LFPmat(:,tt,ch,aS) = tmp(ch,round(trialStarts(tt)-pre)*params.Fs:round(trialStarts(tt)+post)*params.Fs);
            
            
        end
    end
end


%% parsing behavior
trialType = zeros(1,nTrials);
condition = trigs(trigs>=1 & trigs<=27);


%% setting up codes for PSTHs over conflict types.
% These are the correct codes. Double Checked on 20160216
trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)


%% Rasters and PSTHs
for aS2 = 1:length(data)
    % which alignment spot
    switch aS2
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            nTrials = sum(trigs>=1 & trigs<28);
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>=100 & trigs<=105);
            nTrials = 300;
    end
    
    
    for ch = 1:size(data(aS2).channel,2)
        for un = 1:size(data(aS2).channel(ch).unit,2)
            coherenceDataFile = 'COMEUPWITHAFILENAMINGSYSTEM';
            if exist(coherenceDataFile,'file')
                
            else
                display(sprintf('calculating spike-field coherence within trials for\n    LFP channel %s,\n    AP channel %s,\n    unit %d,\n    aligned on %s.'...
                    ,strtrim(LFPlabels{1}),strtrim(unitChanLabels(ch,1:5)),un,alignName));
                
                % caluclate coherence without averaging over trials.
                tic
                [C(:,:,:,un,ch,aS2),phi(:,:,:,un,ch,aS2),~,~,~,t,f] = cohgramcpt(squeeze(LFPmat(:,:,1,aS2)),data(aS2).channel(ch).unit(un).trial, movingWin, params);
                tspec = linspace(-pre,post,length(t));
                toc
            end
            
            % plotting mean cohereograms over trials.
            figure(ch*1000+un)
            for tTy = 1:4
                switch tTy
                    case 1
                        conflictType = 'none';
                    case 2
                        conflictType = 'spatial';
                    case 3
                        conflictType = 'distractor';
                    case 4
                        conflictType = 'both';
                end
                
                % trial-averaged coherogram for current trial type
                Ctt = squeeze(nanmean(C(:,:,trialType(1:nTrials)==tTy,un,ch,aS2),3));
                
                % plotting coheregram for hard trials.
                plotmultipleaxes(tTy,4,1,0.04,ch*1000+un)
                
                % meats
                try
                    imagesc(tspec,f,Ctt');
                catch
                    if isnan(max(max(max(squeeze(C(:,:,:,un,ch,aS2))))))
                        display('coherogram is all NaNs')
                    end
                end
                
                line([0 0],[params.fpass(1) params.fpass(2)],'color','k','linestyle','--','linewidth',2)
                line([mean(rt(trialType(1:nTrials)==tTy)./1000) mean(rt(trialType(1:nTrials)==tTy)./1000)],[params.fpass(1) params.fpass(2)],'color','k','linestyle','-.','linewidth',2)
                
                % deets
                axis xy square tight
                xlim([-0.5 2])
                set(gca, 'linewidth',2,'fontsize',14)
                colormap(jet)
                colorbar
                set(gca, 'linewidth',2,'fontsize',14)
                xlabel('time (seconds)','fontsize',16)
                ylabel('LFP Frequency (Hz)','fontsize',16)
                title(sprintf('%s trials',conflictType))
            end
            
            
            %% TODO: stats, incl shuffle test.
            
            
            %% saving figures.
            figFlag = 'dropbox';
            if strcmp(figFlag,'local')
                if exist(['../../' patientID],'dir')
                    try
                        fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%d_Unit_%d_Coherence_%saligned',patientID,patientID,sessionNum,inclChans(ch),un,alignName);
                        saveas(ch*1000+un,fName, 'pdf')
                        close(ch*1000+un)
                    catch
                        mkdir(sprintf('../../%s/Figs/',patientID))
                        fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%d_Unit_%d_Coherence_%saligned',patientID,patientID,sessionNum,inclChans(ch),un,alignName);
                        saveas(ch*1000+un,fName, 'pdf')
                        close(ch*1000+un)
                    end
                elseif exist('../../Figs','dir')
                    fName = sprintf('../../Figs/%s_session_%d_Channel_%d_Unit_%d_Coherence_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
                    saveas(ch*1000+un,fName, 'pdf')
                    close(ch*1000+un)
                else
                    fName = sprintf('%s_session_%d_Channel_%d_Unit_%d_Coherence_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
                    saveas(ch*1000+un,fName, 'pdf')
                    close(ch*1000+un)
                end
            elseif strcmp(figFlag,'dropbox')
                display('saving to Elliot"s dropbox. Thanks!')
                dbPath = sprintf('/home/elliot/Dropbox/MSITunits_emu/%s/Figs',patientID);
                fName = sprintf('%s/%s_session_%d_Channel_%d_Unit_%d_Coherence_%saligned',dbPath,patientID,sessionNum,inclChans(ch),un,alignName);
                saveas(ch*1000+un,fName, 'pdf')
                close(ch*1000+un)
            end
            
            display(sprintf('figure saved to %s',fName));
            
            %
            %             %% saving stats.
            %             if exist(['./' patientID],'dir')
            %                 try
            %                     fName = sprintf('./%s/Data/%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
            %                     save([fName '.mat'],'neuronStats')
            %                 catch
            %                     mkdir(sprintf('./%s/Data/',patientID))
            %                     fName = sprintf('./%s/Data/%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
            %                     save([fName '.mat'],'neuronStats')
            %                 end
            %             elseif exist('./Data','dir')
            %                 fName = sprintf('./Data/%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
            %                 save([fName '.mat'],'neuronStats')
            %             else
            %                 fName = sprintf('%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
            %                 save([fName '.mat'],'neuronStats')
            %             end
            
            
        end % looping over units for each channel and align spot
    end % looping over channels for each align spot.
end % looping over align spots (Stimulus & response)

