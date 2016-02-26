function [neuronStats] = analyzeMSITunits(patientID,sessionNum,nevFile)
%ANALYZEMSITUNITS plots PSTHs for MSIT units
%
%   [neuronStats] = analyzeMSITunits(patientID,sessionNumber,nevFile) will
%   plot PSTHs for sorted units recorded in a blackrock NEV file or its
%   associated .mat and return statistics.


%% loading data from NEV file
display('loading data...')
% [pathstr, name, ext] = fileparts(nevFile);

ext = nevFile(end-3:end);
if strcmp(ext,'.nev')
    NEV = openNEV(nevFile,'read');
elseif strcmp(ext,'.mat')
    load(nevFile);
end

trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);

%% look for the start of the task
if isequal(trigs(1),255)
    display('task started in this recording')
else
    display('no task start found...')
    display('recording may have started after the behavioral task.')
end


% %% if nev file doesn't contain sorted units:
% if exist(sortedUnits)
%     % getting sorted unit data from offline sorter output (text oer excel).
%     [pathstr, name, ext] = fileparts(sortedUnits);
%     if strcmp(ext,'.txt')
%         Smat = tdfread(sortedUnits)
%         ChanUnitTimestamp = Smat.Channel0x2CUnit0x2CTimestamp0x2CPC_10x2CPC_20x2CPC_3(:,1:3);
%     elseif strcmp(ext,'.xls','xlsx')
%         [Smat,txt,raw] = xlsread(sortedUnits);
%         ChanUnitTimestamp = Smat(:,1:3);
%     end
% end


%% timing (seconds)
pre = 2;
post = 3;


%% creating neural timing variable
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];
ChanUnitTimestamp(ChanUnitTimestamp(:,2)==255,:) = [];

inclChans = unique(ChanUnitTimestamp(:,1));
nChans = length(inclChans);


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
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>100 & trigs<104);
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
                
                %% putting the data in a structure
                data(aS).channel(ch).unit(un).trial(tt).times = unitTimes(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post) - repmat(trialStarts(tt)-pre,length(unitTimes(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post)),1);
                
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


%% parsing behavior
trialType = zeros(1,nTrials);
condition = trigs(trigs>=1 & trigs<=27);


%% setting up codes for PSTHs over conflict types. 
% These are the correct codes. Double Checked on 20160216
trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)


%% setting up codes for button selectivity. 


%% Rasters and PSTHs 
for aSS = 1:length(data)
    % which alignment spot
    switch aSS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>100 & trigs<104);
    end
    
    for ch = 1:size(data(aSS).channel,2)
        for un = 1:size(data(aSS).channel(ch).unit,2)
            
            display('plotting raster over conflict... grab a cocktail, this may take a while.')
            %% plotting rasters and PSTHs
            figure(aSS)
            ah_ras = plotmultipleaxes(1,1,2,0.08,aSS);
            hold on
            for tt = 1:nTrials
                % changing raster color based on trial type
                if trialType(tt)==1
                    rasCol = rgb('DarkGreen');
                elseif trialType(tt)==2
                    rasCol = rgb('Goldenrod');
                elseif trialType(tt)==3
                    rasCol = rgb('OrangeRed');
                elseif trialType(tt)==4
                    rasCol = rgb('FireBrick');
                end
                
                % plotting rasters for conflict (in the least efficient way possible)
                for sp = 1:size(data(aSS).channel(ch).unit(un).trial(tt).times,1)
                    try
                        line([data(aSS).channel(ch).unit(un).trial(tt).times(sp)-pre data(aSS).channel(ch).unit(un).trial(tt).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', rasCol)
                    catch
                        line([data(aSS).channel(ch).unit(un).trial(tt).times(sp)-pre data(aSS).channel(ch).unit(un).trial(tt).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', 'k')
                    end
                end
                
                % stimulus timing lines
                line([0 0], [0 nTrials],'linestyle', '--', 'color', 'k')
                % raster plot details
                xlim([-1 1.5])
                ylim([0 nTrials])
                str = sprintf('patient %s, Channel %d, Unit %d; aligned on %s',patientID ,ch ,un ,alignName);
                title(str,'fontsize',18);
                ylabel('Trials','fontsize', 16)
                set(gca, 'linewidth', 2, 'fontsize', 16);
                
                
                %% plotting rasters for button selectivity. 
                
                
            end
            hold off
            
            
            %% PSTH
            % calculating psths
            kernelWidth = 25  ./1000;
            [Reasy,t,Eeasy] = psth(data(aSS).channel(ch).unit(un).trial(trialType==1), kernelWidth, 'n', [0 pre+post]);
            [Rhard,t,Ehard] = psth(data(aSS).channel(ch).unit(un).trial(trialType==4), kernelWidth, 'n', [0 pre+post]);
            
            
            try
                [Rmedi1,t,Emedi1] = psth(data(aSS).channel(ch).unit(un).trial(trialType==2), kernelWidth, 'n', [0 pre+post]);
            catch
                display('no spatial conflict')
            end
            
            try
                [Rmedi2,t,Emedi2] = psth(data(aSS).channel(ch).unit(un).trial(trialType==3), kernelWidth, 'n', [0 pre+post]);
            catch
                display('no distracter conflict')
            end
            
            %         % overall PSTH
            %         [R,t,E] = psth(data.channel(ch).unit(un).trial, kernelWidth, 'n', [-pre post]);
            tsec = t-repmat(pre,1,length(t));
            
            
            %% generate statistics for each Neuron
            [pEasy,Heasy] = ranksum(Reasy(tsec>=-1 & tsec<=-0.5),Reasy(tsec>=0 & tsec<=1));
            [pHard,Hhard] = ranksum(Rhard(tsec>=-1 & tsec<=-0.5),Rhard(tsec>=0 & tsec<=1));
            neuronStats{ch,un} = table(pEasy,pHard,'VariableNames',{'pEasy','pHard'});
            
            
            %% make it pretty
            % colors and colors and colors and colors
            EcolEasy = rgb('SpringGreen');
            EcolHard = rgb('LightCoral');
            EcolMedi1 = rgb('Orange');
            EcolMedi2 = rgb('Khaki');
            
            % PSTH plot
            ah_psth = plotmultipleaxes(2,1,2,0.08,aSS);
            hold on
            
            
            %% plotting psths
            try
                % plotting med psth
                patch([tsec fliplr(tsec)],[Rmedi1+Emedi1 fliplr(Rmedi1-Emedi1)], EcolMedi1,'edgecolor','none')
                plot(tsec,Rmedi1,'color',rgb('orangered'),'linewidth',2)
            catch
                display('no spatial conflict')
            end
            
            try
                % plotting med psth
                patch([tsec fliplr(tsec)],[Rmedi2+Emedi2 fliplr(Rmedi2-Emedi2)], EcolMedi2,'edgecolor','none')
                plot(tsec,Rmedi2,'color',rgb('goldenrod'),'linewidth',2)
            catch
                display('no spatial conflict')
            end
            
            % plotting easy psth
            patch([tsec fliplr(tsec)],[Reasy+Eeasy fliplr(Reasy-Eeasy)], EcolEasy,'edgecolor','none')
            plot(tsec,Reasy,'color',rgb('DarkGreen'),'linewidth',2)
            
            % plotting hard psth
            patch([tsec fliplr(tsec)],[Rhard+Ehard fliplr(Rhard-Ehard)], EcolHard,'edgecolor','none')
            plot(tsec,Rhard,'color',rgb('darkRed'),'linewidth',2)
            
            % stimulus timing lines
            line([0 0], [0 nTrials],'linestyle', '--', 'color', 'k')
            % PSTH plot details
            xlim([-1 1.5])
            ylim([0 max(Rhard+Ehard)+10])
            hold off
            xlabel('Time (seconds)', 'fontsize', 16);
            ylabel('Firing Rate (spikes/second)', 'fontsize', 16);
            set(gca, 'linewidth', 2, 'fontsize', 16);
            
            
            %% saving figures.
            figFlag = 1;
            if figFlag
                if exist('./Figs','dir')
                    fName = sprintf('./Figs/%s_session_%d_Channel_%d_Unit_%d_Conflict_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
                    saveas(aSS,fName, 'pdf')
                    close(aSS)
                else
                    fName = sprintf('%s_session_%d_Channel_%d_Unit_%d_Conflict_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
                    saveas(aSS,fName, 'pdf')
                    close(aSS)
                end
            end
            
            
            %% plotting PSTHs for button selectivity. 
            
            
            
        end % looping over units for each channel and align spot
    end % looping over channels for each align spot. 
end % looping over align spots (Stimulus & response)

