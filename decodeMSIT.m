function [decodeResults] = decodeMSIT(patientID,sessionNum,nevFile)
%DECODEMSIT classifies trials into conflict conditions based on AP firing.
%
%

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


%% Cues and response times.
cueTimes = trigTimes(trigs>=1 & trigs<=27);
respTimes = trigTimes(trigs>=100 & trigs<=103);

RTs = respTimes-cueTimes;


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


kernel0 = 25  ./1000;
%% Rasters and PSTHs
for aS2 = 1:length(data)
    
    %% initializing decode variables.
    X = cell(nTrials,1);
    classes = cell(nTrials,1);
    
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
    
    % looping over channels, units and trials.
    for ch = 1:size(data(aS2).channel,2)
        for un = 1:size(data(aS2).channel(ch).unit,2)
            display('Calculating PSTHs...')
            [Rfill,t,~] = psth(data(aS2).channel(ch).unit(un).trial, kernel0, 'n', [0 pre+post]);
            for tt = 1:nTrials
                switch trialType(tt)
                    case 1
                        classes{tt} = 'none';
                    case 2
                        classes{tt} = 'spatial';
                    case 3
                        classes{tt} = 'distractor';
                    case 4
                        classes{tt} = 'both';
                end
                
                % calculating psths
                if ~isempty(data(aS2).channel(ch).unit(un).trial(tt).times)
                    [R,~,~] = psth(data(aS2).channel(ch).unit(un).trial(tt), kernel0, 'n', [0 pre+post]);
                else
                    R = zeros(1,length(Rfill));
                end
                
                % time base in seconds.
                tsec = t-repmat(pre,1,length(t));
                
                % stacking data into a matrix.
                X{tt} = cat(2,X{tt},R(tsec>0 & tsec<=RTs(tt)));
                
            end % looping oer trials.
        end % looping over units for each channel and align spot
    end % looping over channels for each align spot.
    
    
    %% reshaping data and dimensionality reduciton
    % finding the length of data over channels and units during each trial.
    for t2 = 1:nTrials
        Xlen(t2) = length(X{t2});
    end
    maxXlen = max(Xlen);
    
    % padding with zeros
    for t3 = 1:nTrials
        padlen = maxXlen-length(X{t3});
        X{t3} = [X{t3} zeros(1,padlen)];
    end
    X = cell2mat(X);
    
    %% dimensionality reduction
    method = 't-SNE';
    display(sprintf('first using %s method to reduce the data dimensionality',method))
    if strcmp(method,'PCA')
        % using PCA
        [coef,~,~] = pca(X');
        clear X
        numCs = 5;
    elseif strcmp(method,'t-SNE')
        no_dims = 2;
        initial_dims = 300;
        perplexity = 30;
        % using t-SNE
        coef = tsne(X, trialType');
        numCs = no_dims;
    else
        coef = X;
    end
    
    % transposing coef for following code.
    if size(coef,1)>size(coef,2)
        coef = coef';
    end
    
    % scatter plot in PCA space.
    figure(1)
    hold on
    Cmap = distinguishable_colors(4);
    for t4 = 1:nTrials
        scatter(coef(1,t4),coef(2,t4),10,Cmap(trialType(t4),:),'filled')
    end
    legend('none','spatial','distractor','both')
    hold off
    
    
    %% Decode using PCs
    %     decodeData = X;
    trainData = coef(1:numCs,1:100)';
    decodeData = coef(1:numCs,101:200)';
    xvalData = coef(1:numCs,201:300)';
    
    
    %% LDA decode for all classes.
    LDAobj = fitcdiscr(trainData,classes(1:100));
    [decodeLabel,~,~] = predict(LDAobj,decodeData);
    xvalLabel = predict(LDAobj,xvalData);
    
    % compare label with classes(101:200)
    decodePerf = sum(strcmpi(decodeLabel,classes(101:200)))./length(decodeLabel);
    xvalPerf = sum(strcmpi(xvalLabel,classes(201:300)))./length(xvalLabel);
    
    
    %% SVM decode over class combinations
    svmTrainClasses = classes(1:100);
    svmDecodeClasses = classes(101:200);
    svmXvalClasses = classes(201:300);
    
    classperms = nchoosek(1:4,2);
    for pms = 1:size(classperms,1)
        % specifying data sets.
        trainIdx = (trialType(1:100)==repmat(classperms(pms,1),1,length(trainData)) | trialType(1:100)==repmat(classperms(pms,2),1,length(trainData)));
        decodeIdx = (trialType(101:200)==repmat(classperms(pms,1),1,length(decodeData)) | trialType(101:200)==repmat(classperms(pms,2),1,length(decodeData)));
        xvalIdx = (trialType(201:300)==repmat(classperms(pms,1),1,length(xvalData)) | trialType(201:300)==repmat(classperms(pms,2),1,length(xvalData)));
        
        % classifying
        SVMobj = fitcsvm(trainData(trainIdx,:),svmTrainClasses(trainIdx),'KernelFunction','gaussian');
        [svmDecodeLabel] = predict(SVMobj,decodeData(decodeIdx,:));
        [svmXvalLabel] = predict(SVMobj,xvalData(xvalIdx,:));
        
        % calculating decode performance.
        svmDecodePerf(pms) = sum(strcmpi(svmDecodeLabel,svmDecodeClasses(decodeIdx)))./length(svmDecodeLabel);
        svmXvalPerf(pms) = sum(strcmpi(svmXvalLabel,svmXvalClasses(xvalIdx)))./length(svmXvalLabel);
    end
    
    
    %% plotting decode performance.
    figure(aS2)
    bar([[decodePerf; xvalPerf].*100 [svmDecodePerf; svmXvalPerf].*100]')
    line([0.5 1.5],[25 25],'color','k','linestyle','--','linewidth',2)
    line([1.5 7.5],[50 50],'color','k','linestyle','--','linewidth',2)
    ylabel('percent correctly classified')
    xlabel('LDA----------------------------------------SVM-------------------------------')
    set(gca,'xticklabel',{'all4','1v2','1v3','1v4','2v3','2v4','3v4'})
    
    
    %% saving figures.
    figFlag = 'local';
    if strcmp(figFlag,'local')
        if exist(['../../' patientID],'dir')
            try
                fName = sprintf('../../%s/Figs/%s_session_%saligned_decodeResults',patientID,patientID,sessionNum,alignName);
                saveas(aS2,fName, 'pdf')
                close(aS2)
            catch
                mkdir(sprintf('../../%s/Figs/',patientID))
                fName = sprintf('../../%s/Figs/%s_session_%saligned_decodeResults',patientID,patientID,sessionNum,alignName);
                saveas(aS2,fName, 'pdf')
                close(aS2)
            end
        elseif exist('./Figs','dir')
            fName = sprintf('./Figs/%s_session_%saligned_decodeResults',patientID,sessionNum,alignName);
            saveas(aS2,fName, 'pdf')
            close(aS2)
        else
            fName = sprintf('%s_session_%saligned_decodeResults',patientID,sessionNum,alignName);
            saveas(aS2,fName, 'pdf')
            close(aS2)
        end
    elseif strcmp(figFlag,'dropbox')
        
    end
    
    %% TODO: output results. 
%     %% saving stats.
%     if exist(['./' patientID],'dir')
%         try
%             fName = sprintf('./%s/Data/%s_session_%d_ConflictStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
%             save([fName '.mat'],'neuronStats')
%         catch
%             mkdir(sprintf('./%s/Data/',patientID))
%             fName = sprintf('./%s/Data/%s_session_%d_ConflictStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
%             save([fName '.mat'],'neuronStats')
%         end
%     elseif exist('./Data','dir')
%         fName = sprintf('./Data/%s_session_%d_ConflictStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
%         save([fName '.mat'],'neuronStats')
%     else
%         fName = sprintf('%s_session_%d_ConflictStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
%         save([fName '.mat'],'neuronStats')
%     end
    
end % looping over align spots (Stimulus & response)

