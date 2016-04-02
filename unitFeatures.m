function fName = unitFeatures(nevFile,patientID,session)
% UNITFEATURES generates summary data for units in a NEV file.
%
%   fName = unitFeatures(nevFile, patientID, session) will generate a
%   pdf that shows unit waveforms and interspike interval histograms for
%   each unit specified in the NEV file.
%
%   The patientID and session arguments are optional.
%
%   The ouput will idicate the location of the pdf file.
%

% author: EHS20160402

close all;

if isequal(nargin,1)
    patientID = [];
    session = 000;
elseif isequal(nargin,2)
    session = 000;
end


%% loading data from NEV file
display('loading unit data...')
% parsing file extension
ext = nevFile(end-3:end);
if strcmp(ext,'.nev')
    NEV = openNEV(nevFile,'read');
elseif strcmp(ext,'.mat')
    load(nevFile);
else % if there isn't an extension specified...
    try
        load([nevFile '.mat']);
    catch
        NEV = openNEV([nevFile '.nev'],'read');
    end
end


%% Looking for task data.
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
timeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);


%% look for the start of the task
if isequal(trigs(1),255)
    taskStartFlag = true;
else
    taskStartFlag = false;
end


%% creating neural timing variable
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./timeRes)'];
notSpikes = ChanUnitTimestamp(:,2)==255;
ChanUnitTimestamp(notSpikes,:) = [];

WFs = NEV.Data.Spikes.Waveform./4;
WFs(:,notSpikes) = [];

WFtime = linspace(1,(size(WFs,1)-1)./timeRes,size(WFs,1))*1000; % waveform timebase in milliseconds.

inclChans = unique(ChanUnitTimestamp(:,1));
nChans = length(inclChans);

channelNames = {NEV.ElectrodesInfo.ElectrodeLabel};

%% analyzing over channels and units.
% looping over Channels
for ch = 1:nChans
    
    % looping over number of units in the AP data
    nUnits = length(unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));
    for un = 1:nUnits
        
        % getting unit times for the current channel and unit.
        currentUnit = ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un;
        unitTimes = ChanUnitTimestamp(currentUnit,3); % in seconds
        unitWFs = WFs(:,currentUnit);
        
        unitWFbar = fliplr(nanmean(unitWFs'));
        unitWFerr = fliplr(std(double(unitWFs')));
        
        %% calucalting average firing rate.
        FRbar = length(unitTimes)/(unitTimes(end)-unitTimes(1));
        
        %% plotting waveforms.
        figure(un)
        plotmultipleaxes(1,2,1,0.1,un)
        hold on
        patch([WFtime fliplr(WFtime)],[unitWFbar+unitWFerr fliplr(unitWFbar-unitWFerr)],[0 0 0],'FaceAlpha',0.3)
        plot(WFtime,unitWFbar,'k','linewidth',2)
        
        set(gca,'fontsize',14,'linewidth',2)
        axis square
        xlabel('time (ms)','fontsize', 16)
        ylabel('mean/std unit waveform (uV)','fontsize', 16)
        
        title(sprintf('channel %s unit %d average firing rate = %d spikes/s',channelNames{ch}(1:5),un,FRbar))
        
        
        %% plotting ISI histogram.
        figure(un)
        plotmultipleaxes(2,2,1,0.1,un)
        histogram(diff(unitTimes))
        
        set(gca,'fontsize',14,'linewidth',2)
        axis square
        xlabel('time ()','fontsize', 16)
        ylabel('AP count','fontsize', 16)
        
        
        %         % loooping over trials
        %         for tt = 1:nTrials
        %
        %             %% putting the data in a structure
        %             data(aS).channel(ch).unit(un).trial(tt).times = unitTimes(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post) - repmat(trialStarts(tt)-pre,length(unitTimes(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post)),1);
        %
        %         end
        
        %% saving stats and figs.
        if ~isempty(patientID)
            if exist(['../../' patientID],'dir')
                try
                    fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%d_Unit_%d_WFs_and_ISIhisto',patientID,patientID,session,inclChans(ch),un)
                    saveas(un,fName, 'pdf')
                    close(un)
                catch
                    mkdir(sprintf('../../%s/Figs/',patientID))
                    fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%d_Unit_%d_WFs_and_ISIhisto',patientID,patientID,session,inclChans(ch),un)
                    saveas(un,fName, 'pdf')
                    close(un)
                end
            elseif exist('../../Figs','dir')
                fName = sprintf('../../Figs/%s_session_%d_Channel_%d_Unit_%d_WFs_and_ISIhisto',patientID,session,inclChans(ch),un)
                saveas(un,fName, 'pdf')
                close(un)
            else
                fName = sprintf('%s_session_%d_Channel_%d_Unit_%d_WFs_and_ISIhisto',patientID,session,inclChans(ch),un)
                saveas(un,fName, 'pdf')
                close(un)
            end
        else
            fName = sprintf('%s_session_%d_Channel_%d_Unit_%d_WFs_and_ISIhisto',patientID,session,inclChans(ch),un)
            saveas(un,fName, 'pdf')
            close(un)
        end
        
    end % looping over units for each channel and align spot
end % looping over channels for each align spot.

