function [behavioralStats,pFeedback,grattonStats] = analyzeMSITbehavior(patientID,sessionNum,nevFile)
% ANALYZEMSITBEHAVIOR analyze bahavior from psychtoolbox MSIT recorded on
%   blackrock.
%
%   [behavioralStats] = analyzeMSITunits(patientID,sessionNumber,nevFile) will
%   plot reaction times over conflict conditions and return statistics.
%
% to fix: previous trial effects & stats

% Author: EHS20160210

% deal with this: ~~~~
saveFlag = 1;

 
%% loading data from NEV file
display('loading data...')
% [pathstr, name, ext] = fileparts(nevFile);
ext = nevFile(end-3:end);
if strcmp(ext,'.nev')
    NEV = openNEV(nevFile);
elseif strcmp(ext,'.mat')
    load(nevFile);
end
trigs = double(NEV.Data.SerialDigitalIO.UnparsedData);
trigTimes = double(NEV.Data.SerialDigitalIO.TimeStampSec);
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);


%% look for the start of the task
if isequal(trigs(1),255)
    display('task started in this recording.')
else
    display('no task start found...')
    display('recording may have started after the behavioral task.')
end


%% parsing behavior
trialType = zeros(1,nTrials);
condition = trigs(trigs>=1 & trigs<28);


%% These are the correct codes. Double Checked on 20160216
trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)


%% Cues and response times. 
cueTimes = trigTimes(trigs>=1 & trigs<=27);
respTimes = trigTimes(trigs>=100 & trigs<=103);

RTs = respTimes-cueTimes;


%% do statistics on RTs over conflict
[P,~,behavioralStats] = anova1(RTs,trialType,'off');

% %sanity check
% sanity1 = RTs(trialType==1)


%% plot RTs over conflict
figure
boxplot(RTs,trialType)
title(['session ' num2str(sessionNum) ', omnibus p-value = ' num2str(P)])
set(gca,'linewidth',2,'fontsize',16,'XTickLabel',{'none','spatial','distractor','both'})
xlabel('Conflict Type')
ylabel('RT (seconds)')
ylim([0 5])
axis square

% saving.
if exist('./Figs','dir')
    savestr = sprintf('./Figs/%s_session%d_conflictBehavior.pdf',patientID,sessionNum)
    saveas(gcf,savestr)
    close(gcf)
else
    savestr = sprintf('%s_session%d_conflictBehavior.pdf',patientID,sessionNum)
    saveas(gcf,savestr)
    close(gcf)
end


%% RTs following feedback.
feedbackMarkers = trigs(trigs>=200 & trigs<=205);
feedbackValence = double(feedbackMarkers<202) + (double(feedbackMarkers>202).*2);
% adjusting for previous trial.
postFeedbackValence =  [NaN; feedbackValence(1:end-1)];


%% do feedBback stats
[pFeedback, ~] = ranksum(RTs(postFeedbackValence==1),RTs(postFeedbackValence==2));

% try
%     RTafterIncorrectFeedback = RTs(feedbackMarkers==201);
% catch
%     display('no incorrect trials')
% end


%% plot RTs following feedback.
figure
boxplot(RTs,postFeedbackValence)
title(['session ' num2str(sessionNum) ', rank sum p-value = ' num2str(pFeedback)])
set(gca,'linewidth',2,'fontsize',16,'XTickLabel',{'Correct','Neutral'})
xlabel('Feedback Valence')
ylabel('RT (seconds)')
ylim([0 5])
axis square

% saving.
if exist('./Figs','dir')
    savestr = sprintf('./Figs/%s_session%d_feedbackBehavior.pdf',patientID,sessionNum)
    saveas(gcf,savestr)
    close(gcf)
else
    savestr = sprintf('%s_session%d_feedbackBehavior.pdf',patientID,sessionNum)
    saveas(gcf,savestr)
    close(gcf)
end


%% RTs for previous trial effects
% conflict = trialType~=1;
% 
% noB4conf = [0 diff(conflict)==1];
% confB4no = [0 diff(conflict)==-1];
% confB4conf = [0 diff(conflict==1)==0];
% noB4no = [[0 diff(conflict==0)==0];
% 
% 
% 
% 
% %% plot RTs for previous trial.
% figure
% boxplot(RTs,gratton)
% title(['session ' num2str(sessionNum) ', omnibus p-value = ' num2str(P)])
% set(gca,'linewidth',2,'fontsize',16,'XTickLabel',{'easy/easy','hard/easy','hard/hard','easy/hard'})
% xlabel('Conflict Type')
% ylabel('RT (seconds)')
% ylim([0 5])
% axis square
% 
% % saving.
% if exist('./Figs','dir')
%     savestr = sprintf('./Figs/%s_session%d_grattonBehavior.pdf',patientID,sessionNum)
%     saveas(gcf,savestr)
%     close(gcf)
% else
%     savestr = sprintf('%s_session%d_grattonBehavior.pdf',patientID,sessionNum)
%     saveas(gcf,savestr)
%     close(gcf)
% end

end

