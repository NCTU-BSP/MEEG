function varargout = process_add_noise_v1( varargin )
% PROCESS_ADD_NOISE: ADD NOISE TO THE INPUT MEASUREMENTS.
%
% USAGE:  OutputFiles = process_simulate_recordings('Run', sProcess, sInputA)
 
% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2014 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Hui-Ling Chan, 2014

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Add noise to recordings';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Simulate'; 
    sProcess.Index       = 916; 
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    % === TIME WINDOW ===
    sProcess.options.label1.Comment = '<HTML><B><U>Input options</U></B>:';
    sProcess.options.label1.Type    = 'label';
    % === Source Active Time Range
    sProcess.options.act_range.Comment = 'Source active range: ';
    sProcess.options.act_range.Type    = 'timewindow';
    sProcess.options.act_range.Value   = [];    
    % === Trials
    sProcess.options.trial_no.Comment = 'Number of trials:';
    sProcess.options.trial_no.Type    = 'value';
    sProcess.options.trial_no.Value   = {1,'',0};     
    
    % background noise
    sProcess.options.label_bn.Comment = '<HTML><BR><B>Background sources</B>';
    sProcess.options.label_bn.Type    = 'label';
    % === number of background sources
    sProcess.options.backsourcenum.Comment = 'Number of sources:';
    sProcess.options.backsourcenum.Type    = 'value';
    sProcess.options.backsourcenum.Value   = {1000,'',0};
    % === STD
    sProcess.options.std.Comment = 'Std of source amplitude:';
    sProcess.options.std.Type    = 'value';
    sProcess.options.std.Value   = {1,'nAm',3};  
    % OPTIONS: orientation type of background sources
    sProcess.options.backsourceori.Comment = {'Normal to cortex','Random','Orientation type:' };
    sProcess.options.backsourceori.Type    = 'radio_line';
    sProcess.options.backsourceori.Value   = 2;
  
    
    % sensor noise
    sProcess.options.label_sn.Comment = '<HTML><BR><B>Sensor noise</B>';
    sProcess.options.label_sn.Type    = 'label';
    % OPTIONS: orientation type of background sources
    sProcess.options.sensornoisetype.Comment = {'Fixed weights','Fixed SNR','None','Strength:'};
    sProcess.options.sensornoisetype.Type    = 'radio_line';
    sProcess.options.sensornoisetype.Value   = 1;
    % === weight
    sProcess.options.sensornoiseweight.Comment = 'Fixed weights: ';
    sProcess.options.sensornoiseweight.Type    = 'value';
    sProcess.options.sensornoiseweight.Value   = {1,'',3};   
    % === SNR
    sProcess.options.SNR.Comment = 'Fixed SNR: ';
    sProcess.options.SNR.Type    = 'value';
    sProcess.options.SNR.Value   = {1,'dB',1};    


%     % === std of background sources
%     sProcess.options.backsourcestd.Comment = 'Std:';
%     sProcess.options.backsourcestd.Type    = 'value';
%     sProcess.options.backsourcestd.Value   = {1,'dB',2};

    % sensor noise
    sProcess.options.label_rw.Comment = '<HTML><BR>';
    sProcess.options.label_rw.Type    = 'label';
    sProcess.options.savenoise.Comment = 'Save noise recordings';
    sProcess.options.savenoise.Type    = 'checkbox';
    sProcess.options.savenoise.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    OutputFiles = {}; 
    nTrials = length(sInput);
    bst_progress('start', 'Applying process: Add noise to recordings', '', 0, nTrials*3);
    %% =================
    % === GET OPTIONS ===
    BSnum = sProcess.options.backsourcenum.Value{1};   
    BSstd = sProcess.options.std.Value{1};
    BSori = sProcess.options.backsourceori.Value;
    SNtype = sProcess.options.sensornoisetype.Value;
    SNweight = sProcess.options.sensornoiseweight.Value{1};
    SNR = sProcess.options.SNR.Value{1};
    isSaveNoise = sProcess.options.savenoise.Value;
    SourceActiveTimeRange = sProcess.options.act_range.Value{1};
    TrialNo = sProcess.options.trial_no.Value{1};
    
    InputsData = sInput ;
    % ===== LOAD CHANNEL FILE =====
    % Get condition
    sStudy = bst_get('Study', sInput.iStudy);
    % Load channel file
    ChannelMat = in_bst_channel(InputsData(1).ChannelFile);
    
    % ===== LOAD HEAD MODEL =====
    % Get channel study
    [sChannelStudy, iChannelStudy] = bst_get('ChannelFile', InputsData(1).ChannelFile);
    % Load the default head model
    HeadModelFile = sChannelStudy.HeadModel(sChannelStudy.iHeadModel).FileName;
    sHeadModel = in_headmodel_bst(HeadModelFile);
    
    % Get all the MEG/EEG channels
    Modalities = {};
    if ~isempty(sHeadModel.MEGMethod)
        Modalities{end+1} = 'MEG';
    end
    if ~isempty(sHeadModel.EEGMethod)
        Modalities{end+1} = 'EEG';
    end
    if ~isempty(sHeadModel.SEEGMethod)
        Modalities{end+1} = 'SEEG';
    end
    if ~isempty(sHeadModel.ECOGMethod)
        Modalities{end+1} = 'ECOG';
    end
    iChannels = channel_find(ChannelMat.Channel, Modalities);
       
    % Get forward field
    Kernel = sHeadModel.Gain(iChannels,:);
    %Kernel(abs(Kernel(:)) < eps) = eps; % Set zero elements to strictly non-zero
    [nChannels ,nSources] = size(Kernel); % size of Gain Matrix
    
    % ===== LOAD NOISE COVARIANCE =====
    if SNtype == 1 ||  SNtype == 2
        NoiseCovMat  = load(file_fullpath(sChannelStudy.NoiseCov.FileName));
        if isempty(NoiseCovMat.NoiseCov)
            bst_report('Error', sProcess, [], 'No noise covariance matrix available in this study.');
            return;
        end
        NoiseCov     = NoiseCovMat.NoiseCov(iChannels,iChannels);
    end
       
    
    for i = 1:nTrials
        
        DataMat = in_bst(InputsData(i).FileName, [], 0);
        iSourceActiveTime = panel_time('GetTimeIndices', DataMat.Time, SourceActiveTimeRange);
        dataP = TrialNo*mean(mean(DataMat.F(iChannels,iSourceActiveTime).^2));
        nTime = length(DataMat.Time);
        % === GET BACKGROUND SOURCES ===
        % Deploy random sources for generating background activity
        if BSnum > 0
            bSourceGridIndex = randi(nSources/3,[1 BSnum]);
            bSourceGridIndex3 = reshape([(bSourceGridIndex-1)*3+1; (bSourceGridIndex-1)*3+2; (bSourceGridIndex-1)*3+3],BSnum*3,1);
            % bSourceSignal = wgn(nBSources,size(Sources.ImageGridAmp,2),powerBSources)*(10^-12);
            bSourceSignal = randn(BSnum,nTime)*BSstd*(1e-9);

            if BSori == 2
                % Generate random orientation for background sources 
                bSourceOri = zeros(BSnum,3);
                bSourceOri(:,3) = rand(BSnum, 1)*2-1;
                t = rand(BSnum, 1)*2*pi;
                r = sqrt(1-(bSourceOri(:,3).^2));
                bSourceOri(:,1) = r.*cos(t);
                bSourceOri(:,2) = r.*sin(t);
                bSourceOri = bSourceOri./norm(bSourceOri);
            else
                % If no orientations: error
                if isempty(sHeadModel.GridOrient)
                    bst_report('Error', sProcess, [], 'No source orientations available in this head model.');
                    return;
                end
                bSourceOri = sHeadModel.GridOrient(bSourceGridIndex,:);
            end

            % Calculate Lead field vector for each background source 
            bSourceLeadfield = bst_gain_orient(sHeadModel.Gain(iChannels,bSourceGridIndex3), bSourceOri);

            % Calculate measurements contributed by background sources 
            bF = bSourceLeadfield * bSourceSignal;
        else
            bF = zeros(length(iChannels),nTime);
        end
        bst_progress('inc',1);
        % === GET SENSOR NOISE ===
        if SNtype == 1 ||  SNtype == 2
            %R = chol(noise.NoiseCov(Sources.GoodChannel,Sources.GoodChannel));
            R = (diag(NoiseCov)).^0.5;
            cF = repmat(R,1,nTime).*randn(length(iChannels),nTime);  
            
            if SNtype == 2
                mxRatio = 10;
                mnRatio = 0;              
                generatedSNR = 100;
                lastGeneratedSNR = 0;
                while abs(generatedSNR - SNR) > abs(SNR*0.05) 
                    ra = (mxRatio + mnRatio)/2;
                    noiseF = bF + ra*cF;
                    generatedSNR = 10*log10(dataP/mean(noiseF(:).^2));
                    if generatedSNR > SNR*1.05                       
                        mnRatio = ra;
                    elseif generatedSNR < SNR*0.95
                        mxRatio = ra;
                    end 
%                     disp(['SigGen> Noise level ' num2str(generatedSNR) 'dB (ratio = ' num2str(ra) ')']);
                    if abs(lastGeneratedSNR-generatedSNR) < abs(SNR*0.01) 
                        bst_report('Error', sProcess, [], 'The amplitude of background noise is too large or small.');
                        return;
                    end
                    lastGeneratedSNR = generatedSNR;
                end
            elseif SNtype == 1
                noiseF = bF + SNweight*cF;
            end
        else
            %cF = zeros(length(Sources.GoodChannel),size(Sources.ImageGridAmp,2));
            noiseF = bF;        
        end

        sF = DataMat.F(iChannels,:);
        finalSNR = 10*log10(dataP./[TrialNo 1]/mean(noiseF(:).^2));
        disp(['SigGen> Noise level ' num2str(finalSNR(1)) 'dB (' num2str(TrialNo) ' trials:' num2str(finalSNR(2)) 'dB)']);
        
        bst_progress('inc',1);
        % === SAVE RECORDINGS ===
        % Generate data matrix
        F = DataMat.F;
        F(iChannels,:) = sF + noiseF;

        % Create a new data file structure
        outDataMat = db_template('datamat');
        outDataMat.F           = F;
        outDataMat.Comment     = [DataMat.Comment '| add noise'];
        outDataMat.ChannelFlag = ones(length(ChannelMat.Channel), 1);
        outDataMat.Time        = DataMat.Time;
        outDataMat.DataType    = 'recordings';
        outDataMat.Device      = 'simulation';
        outDataMat.nAvg        = 1;
        outDataMat.Events      = [];
        outDataMat.SNR         = finalSNR;
        % Add history entry
        outDataMat = bst_history('add', outDataMat, 'simulate', ['Add noise to file: ' InputsData(i).FileName]);
        % Output filename
        DataFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sim');
        % Save on disk
        bst_save(DataFile, outDataMat, 'v6');
        % Register in database
        db_add_data(sInput.iStudy, DataFile, outDataMat);
        % Return data file
        OutputFiles{i} = DataFile; 
        bst_progress('inc',1);
        
        if isSaveNoise
            % === SAVE RECORDINGS ===
            % Generate data matrix
            F = DataMat.F;
            F(iChannels,:) = noiseF;

            % Create a new data file structure
            outDataMat = db_template('datamat');
            outDataMat.F           = F;
            outDataMat.Comment     = ['noise' ];
            outDataMat.ChannelFlag = ones(length(ChannelMat.Channel), 1);
            outDataMat.Time        = DataMat.Time;
            outDataMat.DataType    = 'recordings';
            outDataMat.Device      = 'simulation';
            outDataMat.nAvg        = 1;
            outDataMat.Events      = [];
            % Add history entry
            outDataMat = bst_history('add', outDataMat, 'simulate', ['Add noise to file: ' InputsData(i).FileName]);
            % Output filename
            DataFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sim');
            % Save on disk
            bst_save(DataFile, outDataMat, 'v6');
            % Register in database
            db_add_data(sInput.iStudy, DataFile, outDataMat);
            % Return data file
            OutputFiles{i} = DataFile;
        end
        bst_progress('inc',1);
        clear DataMat;
        clear F;
        clear outDataMat;
    end
    bst_progress('stop');
end
