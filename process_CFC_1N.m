function varargout = process_pac1N( varargin )
% PROCESS_BEAMFORMER_TEST: 

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
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
% Authors: 

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Phase-amplitude coupling 1xN - Onslow';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Connectivity';
    sProcess.Index       = 760;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw',      'data',     'results',  'matrix'};
    sProcess.OutputTypes = {'timefreq', 'timefreq', 'timefreq', 'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    % Separator
    %sProcess.options.sep.Type     = 'separator';
    sProcess.options.sep.Comment = '  ';
    sProcess.options.sep.Type    = 'label';

    sProcess.options.result_comm.Comment    = 'Comment: ';
    sProcess.options.result_comm.Type       = 'text';
    sProcess.options.result_comm.Value      = '';
%     
%     sProcess.options.ref_lag.Comment = 'Time lag (+, reference leading): ';
%     sProcess.options.ref_lag.Type    = 'value';
%     sProcess.options.ref_lag.Value   = {0, 'ms', 1};
% 
%     % === ACTIVE TIME RANGE
%     sProcess.options.corr_range.Comment = 'Time range of interest: ';
%     sProcess.options.corr_range.Type    = 'timewindow';
%     sProcess.options.corr_range.Value   = [];
    % === ACTIVE TIME WINDOW SIZE
%     sProcess.options.active_window_size.Comment = 'Sliding window size: ';
%     sProcess.options.active_window_size.Type    = 'value';
%     sProcess.options.active_window_size.Value   = {0.02, 'ms', 1};
%     % === ACTIVE TEMPORAL RESOLUTION
%     sProcess.options.corr_tresolution.Comment = 'Temporal resolution: ';
%     sProcess.options.corr_tresolution.Type    = 'value';
%     sProcess.options.corr_tresolution.Value   = {0.01, 'ms', 1};

    % ==== ESTIMATOR ====
    sProcess.options.label_pac.Comment = '<HTML><BR><B><U>Estimator options</U></B>:';
    sProcess.options.label_pac.Type    = 'label';
    % === PAC MEASURE ===
    sProcess.options.pacmeasure.Comment    = {'ESC', 'MI', 'CFC', 'Pac measure:'};
    sProcess.options.pacmeasure.Type       = 'radio_line';
    sProcess.options.pacmeasure.Value      = 1;
    % === TF METHOD  ===
    sProcess.options.tfmethod.Comment    = {'Hilbert', 'Wavelet', 'STFT', 'TF method:'};
    sProcess.options.tfmethod.Type       = 'radio_line';
    sProcess.options.tfmethod.Value      = 1;
    % === WINDOW LENGTH
%     sProcess.options.winlength.Comment = 'Estimator window length: ';
%     sProcess.options.winlength.Type    = 'value';
%     sProcess.options.winlength.Value   = {0.128, 'ms ', 1};
%     sProcess.options.winlength.InputTypes = {'data'};
    % === Overlap
    sProcess.options.winoverlap.Comment = 'Overlap percentage: ';
    sProcess.options.winoverlap.Type    = 'value';
    sProcess.options.winoverlap.Value   = {0.75, '% ', 1};
%     sProcess.options.winoverlap.InputTypes = {'data'};
    % Options: Time-freq
    % === NESTING FREQ
    sProcess.options.nesting.Comment = 'Nesting frequency band (low):';
    sProcess.options.nesting.Type    = 'range';
    sProcess.options.nesting.Value   = {[2, 30], 'Hz', 2};
    % === FREQ 
    sProcess.options.nestingwidth.Comment = 'Nesting frequency step (low):';
    sProcess.options.nestingwidth.Type    = 'value';
    sProcess.options.nestingwidth.Value   = {0.75, 'Hz', 2};  
    % === NESTED FREQ
    sProcess.options.nested.Comment = 'Nested frequency band (high):';
    sProcess.options.nested.Type    = 'range';
    sProcess.options.nested.Value   = {[40, 150], 'Hz', 2};
    % === FREQ 
    sProcess.options.width.Comment = 'Cycles of wavelet:';
    sProcess.options.width.Type    = 'value';
    sProcess.options.width.Value   = {7, 'Cycles', 0};   
    % === Freq band
    sProcess.options.freqband.Comment = 'Frequency band of interested:';
    sProcess.options.freqband.Type    = 'text';
    sProcess.options.freqband.Value   = 'beta';
    sProcess.options.freqband.InputTypes = {'timefreq'};
    % === TIME LAGGED
    sProcess.options.tlag.Comment = 'Time-lagged';
    sProcess.options.tlag.Type    = 'checkbox';
    sProcess.options.tlag.Value   = 1;
    sProcess = DefineConnectOptions(sProcess, 0);
    
    % ==== OUTPUT ====
    sProcess.options.label_out.Comment = '<HTML><BR><U><B>Output configuration</B></U>:';
    sProcess.options.label_out.Type    = 'label';
    % === AVERAGE OUTPUT FILES
    sProcess.options.avgoutput.Comment = 'Save average PAC across trials (one output file only)';
    sProcess.options.avgoutput.Type    = 'checkbox';
    sProcess.options.avgoutput.Value   = 1;
    
    % === CONNECT INPUT
    %sProcess = process_corr1n('DefineConnectOptions', sProcess, 1);
end

%% ===== DEFINE SCOUT OPTIONS =====
function sProcess = DefineConnectOptions(sProcess, isConnNN) %#ok<DEFNU>
    % === TIME WINDOW ===
    sProcess.options.label1.Comment = '<HTML><B><U>Input options</U></B>:';
    sProcess.options.label1.Type    = 'label';
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % === FROM: CONNECTIVITY [1xN] ===
    if ~isConnNN
        % === FROM: REFERENCE CHANNELS ===
        sProcess.options.src_channel.Comment    = 'Source channel: ';
        sProcess.options.src_channel.Type       = 'channelname';
        sProcess.options.src_channel.Value      = 'name';
        sProcess.options.src_channel.InputTypes = {'data','raw'};
        % === FROM: ROW NAME ===
        sProcess.options.src_rowname.Comment    = 'Source rows (names or indices): ';
        sProcess.options.src_rowname.Type       = 'text';
        sProcess.options.src_rowname.Value      = '';
        sProcess.options.src_rowname.InputTypes = {'timefreq', 'matrix'};
    end
    % === TO: SENSOR SELECTION ===
    sProcess.options.dest_sensors.Comment    = 'Sensor types or names (empty=all): ';
    sProcess.options.dest_sensors.Type       = 'text';
    sProcess.options.dest_sensors.Value      = 'MEG, EEG';
    sProcess.options.dest_sensors.InputTypes = {'data','raw'};
    % === SCOUTS ===
    sProcess.options.scouts.Comment = 'Use scouts';
    if isConnNN
        sProcess.options.scouts.Type = 'scout_confirm';
    else
        sProcess.options.scouts.Type = 'scout';
    end
    sProcess.options.scouts.Value      = [];
    sProcess.options.scouts.InputTypes = {'results'};
    % Atlas: surface/volume
    sProcess.options.isvolume.Comment = '';
    sProcess.options.isvolume.Type    = 'checkbox';
    sProcess.options.isvolume.Value   = 0;
    sProcess.options.isvolume.Hidden  = 1;
    % === SCOUT FUNCTION ===
    sProcess.options.scoutfunc.Comment    = {'Mean', 'Max', 'PCA', 'Std', 'All', 'Scout function:'};
    sProcess.options.scoutfunc.Type       = 'radio_line';
    sProcess.options.scoutfunc.Value      = 2;
    sProcess.options.scoutfunc.InputTypes = {'results'};
    % === SCOUT TIME ===
    sProcess.options.scouttime.Comment    = {'Before', 'After', 'When to apply the scout function:'};
    sProcess.options.scouttime.Type       = 'radio_line';
    sProcess.options.scouttime.Value      = 2;
    sProcess.options.scouttime.InputTypes = {'results'};
end
%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputA) %#ok<DEFNU>

    %% ===== DEFAULT OPTIONS =====
    Def_OPTIONS.Method        = 'corr';
    Def_OPTIONS.ProcessName   = '';
    Def_OPTIONS.TargetA       = [];
    Def_OPTIONS.TargetB       = [];
    Def_OPTIONS.Freqs         = 0;
    Def_OPTIONS.TimeWindow    = [];
    Def_OPTIONS.ScoutFunc     = 'all';         % Scout function {mean, max, pca, std, all}
    Def_OPTIONS.ScoutTime     = 'before';      % When to apply scout function: {before, after}
    Def_OPTIONS.RemoveMean    = 1;             % Option for Correlation
    Def_OPTIONS.CohMeasure    = 'mscohere';    % {'mscohere'=Magnitude-square, 'icohere'=Imaginary}
    Def_OPTIONS.MaxFreqRes    = [];            % Option for spectral estimates (Coherence, spectral Granger)
    Def_OPTIONS.MaxFreq       = [];            % Option for spectral estimates (Coherence, spectral Granger)
    Def_OPTIONS.CohOverlap    = 0.50;          % Option for Coherence
    Def_OPTIONS.GrangerOrder  = 10;            % Option for Granger causality
    Def_OPTIONS.GrangerDir    = 'out';         % Option for Granger causality
    Def_OPTIONS.RemoveEvoked  = 0;             % Removed evoked response to each single trial (useful to bring signals closer to a stationnary state)
    Def_OPTIONS.isMirror      = 1;             % Option for PLV
    Def_OPTIONS.isSymmetric   = 0;            % Optimize processing and storage for simple matrices
    Def_OPTIONS.pThresh       = 0.05;          % Significativity threshold for the metric
    Def_OPTIONS.OutputMode    = 'input';       % {'avg','input','concat'}
    Def_OPTIONS.iOutputStudy  = [];

    % Copy default options to OPTIONS structure (do not replace defined values)
    OPTIONS = struct_copy_fields(sProcess.options, Def_OPTIONS, 0);

    % Initialize returned list of files
    OutputFiles = {};
    OPTIONS.isSymmetric   = 0;
    OPTIONS.ProcessName   = 'PAC';
    % ===== GET OPTIONS =====
    if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
        OPTIONS.TimeWindow = sProcess.options.timewindow.Value{1};
    else
        OPTIONS.TimeWindow = [];
    end
    % Get and check frequencies
    OPTIONS.BandNesting = sProcess.options.nesting.Value{1};
    OPTIONS.BandNested  = sProcess.options.nested.Value{1};
    OPTIONS.Width = sProcess.options.width.Value{1};
    
    if (min(OPTIONS.BandNesting) < 0.5)
        bst_report('Error', sProcess, [], 'This function cannot be used to estimate PAC for nesting frequencies below 1Hz.');
        return;
    end
    if (max(OPTIONS.BandNesting) > min(OPTIONS.BandNested))
        bst_report('Error', sProcess, [], 'The low and high frequency band cannot overlap.');
        return;
    end
    FreqSteps = sProcess.options.nestingwidth.Value{1};
    if FreqSteps == 0
        OPTIONS.Freqs = mean(OPTIONS.BandNesting);
        FreqSteps = OPTIONS.BandNesting(2) - OPTIONS.BandNesting(1);
    else
        OPTIONS.Freqs = (OPTIONS.BandNesting(1)):FreqSteps:(OPTIONS.BandNesting(2));
    end
    % Get target
    if ismember(sInputA(1).FileType, {'data','raw'}) && isfield(sProcess.options, 'src_channel') && ~isempty(sProcess.options.src_channel.Value)
        OPTIONS.Target = sProcess.options.src_channel.Value;
    elseif strcmpi(sInputA(1).FileType, 'results') && isfield(sProcess.options, 'scouts') && ~isempty(sProcess.options.scouts.Value)
        OPTIONS.Target = sProcess.options.scouts.Value;
    elseif ismember(sInputA(1).FileType, {'timefreq', 'matrix'}) && isfield(sProcess.options, 'src_rowname') && ~isempty(sProcess.options.src_rowname.Value)
        OPTIONS.Target = sProcess.options.src_rowname.Value;
    else
        OPTIONS.Target = [];
    end 
    OPTIONS.isTimeLag  = sProcess.options.tlag.Value;
    switch (sProcess.options.pacmeasure.Value)
        case 1, OPTIONS.Method = 'esc';
        case 2, OPTIONS.Method = 'mi';
        case 3, OPTIONS.Method = 'cfc';
    end
    switch (sProcess.options.tfmethod.Value)
        case 1, OPTIONS.TFmethod = 'hilbert';
        case 2, OPTIONS.TFmethod = 'wavelet';
        case 3, OPTIONS.TFmethod = 'stft';
    end
    measure = OPTIONS.Method ;
    OPTIONS.isAvgOutput = sProcess.options.avgoutput.Value;
    if (length(sInputA) == 1)
        OPTIONS.isAvgOutput = 0;
    end
    
    if strcmpi(measure,'cfc')
        SegmentOverlap = sProcess.options.winoverlap.Value{1}/100;
        %SegmentLength = sProcess.options.winlength.Value{1};
    end
    
    % ===== GET SCOUTS OPTIONS =====
    if strcmpi(sInputA(1).FileType, 'results') && isfield(sProcess.options, 'scouts') && isfield(sProcess.options.scouts, 'Value')
        % Selected scouts
        sScouts = sProcess.options.scouts.Value;
        % Override scouts function
        switch (sProcess.options.scoutfunc.Value)
            case 1, OPTIONS.ScoutFunc = 'mean';
            case 2, OPTIONS.ScoutFunc = 'max';
            case 3, OPTIONS.ScoutFunc = 'pca';
            case 4, OPTIONS.ScoutFunc = 'std';
            case 5, OPTIONS.ScoutFunc = 'all';
        end
        % Scout function order
        switch (sProcess.options.scouttime.Value)
            case 1, OPTIONS.ScoutTime = 'before';
            case 2, OPTIONS.ScoutTime = 'after';
        end
        % Perform some checks
        if strcmpi(OPTIONS.ScoutTime, 'before') && ismember(OPTIONS.ScoutFunc, {'max', 'std'})
           bst_report('Error', sProcess, [], 'Scout functions MAX and STD should not be applied before estimating the PAC.');
           return;
        end
        if strcmpi(OPTIONS.ScoutTime, 'after') && strcmpi(OPTIONS.ScoutFunc, 'pca')
            bst_report('Error', sProcess, [], 'Scout function PCA cannot be applied after estimating the PAC.');
            return;
        end
        % Set input/output scouts functions
        if ~isempty(sScouts)
            OPTIONS.Target = sScouts;
            OPTIONS.isScout = 1;
            % Apply function before: get all the scouts time series in advance
            if strcmpi(OPTIONS.ScoutTime, 'before')
                [OPTIONS.TargetFunc] = deal(OPTIONS.ScoutFunc);
            % Apply function after: Get all the time series of all the scouts
            elseif strcmpi(OPTIONS.ScoutTime, 'after')
                [OPTIONS.TargetFunc] = deal('all');
            end
        end
    else
        OPTIONS.isScout = 0;
    end
    
    nAvg = 0;
    % Initialize progress bar
    if bst_progress('isVisible')
        startValue = bst_progress('get');
    else
        startValue = 0;
    end
    % Options for LoadInputFile()
    if strcmpi(sInputA(1).FileType, 'results')
        LoadOptions.LoadFull = 0;  % Load kernel-based results as kernel+data
    else
        LoadOptions.LoadFull = 1;  % Load the full file
    end
    LoadOptions.IgnoreBad   = 1;  % From raw files: ignore the bad segments
    LoadOptions.ProcessName = func2str(sProcess.Function);
    
    
    maxlag = 0;
    DirectPAC_avg = [];
    for iFile = 1:length(sInputA)  
        DirectPAC = [];
        
        % ===== LOAD SIGNALS =====
        bst_progress('text', sprintf('PAC: Loading input file (%d/%d)...', iFile, length(sInputA)));
        
        % Load input signals 
        [sInputRef, nSignalsRef, iRowsRef] = bst_process('LoadInputFile', sInputA(iFile).FileName, OPTIONS.Target, OPTIONS.TimeWindow, LoadOptions);
        if isempty(sInputRef) || isempty(sInputRef.Data)
            return;
        end
        [sInput, nSignals, iRows] = bst_process('LoadInputFile', sInputA(iFile).FileName, [], OPTIONS.TimeWindow, LoadOptions);
        if isempty(sInput) || isempty(sInput.Data)
            return;
        end
        
        % Get sampling frequency
        sRate = 1 / (sInput.Time(2) - sInput.Time(1));
        % Check the nested frequencies
        if (OPTIONS.BandNested(2) > sRate/3)
            % Warning
            strMsg = sprintf('Higher nesting frequency is too high (%d Hz) compared with sampling frequency (%d Hz): Limiting to %d Hz', round(OPTIONS.BandNested(2)), round(sRate), round(sRate/3));
            disp([10 'process_pac> ' strMsg]);
            bst_report('Warning', 'process_pac', [], strMsg);
            % Fix higher frequencyy
            OPTIONS.BandNested(2) = sRate/3;
        end
        % Check the extent of bandNested band
        if (OPTIONS.BandNested(2) <= OPTIONS.BandNested(1))
            bst_report('Error', 'process_pac', [], sprintf('Invalid frequency range: %d-%d Hz', round(OPTIONS.BandNested(1)), round(OPTIONS.BandNested(2))),'n');
            continue;
        end
        
        if ~isempty(sInput.ImagingKernel)
            Fblock = sInput.ImagingKernel * sInput.Data;
        else
            Fblock = sInput.Data;            
        end
        
        if ~isempty(sInputRef.ImagingKernel)            
            FblockRef = sInputRef.ImagingKernel * sInputRef.Data;
        else           
            FblockRef = sInputRef.Data;
        end
        nT = size(Fblock,2);
        nRef = nSignalsRef;
        nComponents = 1;
        if strcmpi(sInput.DataType, 'results')
            nComponents = sInput.nComponents;           
            if nComponents == 0
                error('PAC metrics are not supported for mixed source models.');
            end
            nRef = nSignalsRef/nComponents;
        end
        
        nFreqs = length(OPTIONS.Freqs);
        for iSigX = 1:nComponents:nSignalsRef

            sigX = amp_vec(FblockRef(iSigX+(0:nComponents-1),:),OPTIONS.BandNested,sRate,OPTIONS.Width,OPTIONS.TFmethod);    
            %sigX = ampvec(sum(OPTIONS.BandNested)/2,FblockRef(iSigX+(0:nComponents-1),:),sRate,OPTIONS.BandNested(2)-OPTIONS.BandNested(1));
            if strcmp(measure, 'cfc')
                sigX = (sigX / sRate).^2;
            end
            for iFreq = 1:nFreqs
                if OPTIONS.isTimeLag == 1                        
                    maxlag = 100;%min(floor(sRate/OPTIONS.Freqs(iFreq)/2),floor(nT*0.8)/4);
                end
                BandBounds = [OPTIONS.Freqs(iFreq)-(FreqSteps/2) OPTIONS.Freqs(iFreq)+(FreqSteps/2)];
                nS = 0;
                
                if strcmp(measure, 'esc')
                    sigYmat = bp_vec(Fblock,BandBounds,sRate,OPTIONS.Width,OPTIONS.TFmethod);
                elseif strcmp(measure, 'mi')
                    sigYmat = ph_vec(Fblock,BandBounds,sRate,OPTIONS.Width,OPTIONS.TFmethod);
                end


                for iSigY = 1:nComponents:nSignals
                    
                    if mod(iSigY,round(nSignals/nComponents/100))==0
                        bst_progress('set', round(startValue + (((iFile-1)*nRef*nFreqs+(iSigX-1)*nFreqs+(iFreq-1))*10+nS) / length(sInputA) / nRef /nFreqs /10* 100)); 
                        nS = nS + 1;
                    end
                    
                    %pacmat = zeros(nComponents,nComponents);
                    
                    if strcmp(measure, 'esc')
                        %sigY = bp_vec(Fblock(iSigY+(0:nComponents-1),:),BandBounds,sRate,OPTIONS.Width,OPTIONS.TFmethod);
                        %sigY = bpvec(sum(BandBounds)/2,Fblock(iSigY+(0:nComponents-1),:),sRate,BandBounds(2)-BandBounds(1));
                        %OPTIONS.RemoveMean = 1;
                        %pacmat = bst_corrn(sigX, sigY, OPTIONS.RemoveMean);
                        pacmat = lagged_corr(sigX,sigYmat(iSigY+(0:nComponents-1),:),maxlag);

                    elseif strcmp(measure, 'mi')
                        %sigY = ph_vec(Fblock(iSigY+(0:nComponents-1),:),BandBounds,sRate,OPTIONS.Width,OPTIONS.TFmethod);
                        %sigY = phasevec(sum(BandBounds)/2,Fblock(iSigY+(0:nComponents-1),:),sRate,BandBounds(2)-BandBounds(1));
                        pacmat = mi_measure(sigYmat(iSigY+(0:nComponents-1),:),sigX);
                    elseif strcmp(measure, 'cfc')
                        sigY = Fblock(iSigY+(0:nComponents-1),:);                     
                        pacmat = coherence(sigX, sigY, BandBounds, sRate, OPTIONS.Width, SegmentOverlap);
                    end
                    

                    if isempty(DirectPAC)
                        sInputRef.RowNames = sInputRef.RowNames(1:nComponents:nSignalsRef);
                        sInput.RowNames = sInput.RowNames(1:nComponents:nSignals);
                        DirectPAC =  zeros(nRef, nSignals/nComponents,1, nFreqs);                        
                    end

%                     if nComponents > 1
%                         [u,s,v] = svd(real(pacmat));
%                         csd = s(1,1);
%                     end

                    DirectPAC(ceil(iSigX/nComponents),ceil(iSigY/nComponents),1,iFreq)=pacmat;
                end
            end
            
            
        end
        
        % ===== PROCESS SCOUTS =====
        % If the scout function has to be applied AFTER the PAC computation
        if ~isempty(OPTIONS.Target) && isstruct(OPTIONS.Target) && strcmpi(OPTIONS.ScoutTime, 'after') && ~strcmpi(OPTIONS.ScoutFunc, 'all')
            nScouts = length(OPTIONS.Target);
            DirectPAC_scouts = zeros(nScouts, size(DirectPAC,2), size(DirectPAC,3), size(DirectPAC,4));
            iVerticesAll = [1, cumsum(cellfun(@length, {OPTIONS.Target.Vertices})) + 1];
            % For each unique row name: compute a measure over the clusters values
            for iScout = 1:nScouts
                iScoutVert = iVerticesAll(iScout):iVerticesAll(iScout+1)-1;
                F = reshape(DirectPAC(iScoutVert,:,:,:), length(iScoutVert), []);
                F = bst_scout_value(F, OPTIONS.ScoutFunc);
                DirectPAC_scouts(iScout,:,:,:) = reshape(F, [1, size(DirectPAC,2), size(DirectPAC,3), size(DirectPAC,4)]);
            end
            % Save only the requested rows
            sInput.RowNames = {OPTIONS.Target.Label};
            DirectPAC = DirectPAC_scouts;
        end
        DirectPAC = reshape(DirectPAC, [], 1, nFreqs);
       
        result_comment = sProcess.options.result_comm.Value;
        if ~isempty(result_comment)
            result_comment = [result_comment ':'];
        end
        % Base comment
        if  OPTIONS.isTimeLag
            Comment = [ result_comment 'Lagged' upper(measure) '(' upper(OPTIONS.TFmethod)  ',1xN)'];
        else
            Comment = [ result_comment upper(measure) '(' upper(OPTIONS.TFmethod)  ',1xN)'];
        end

        
        % Time window (RAW only)
        if ~isempty(strfind(sInputA(iFile).Condition, '@raw'))
            Comment = [Comment, sprintf('(%ds-%ds)', round(OPTIONS.TimeWindow))];
        end
        % Scouts
        if isstruct(OPTIONS.Target) && (length(OPTIONS.Target) < 6)
            Comment = [Comment, ':'];
            for is = 1:length(OPTIONS.Target)
                Comment = [Comment, ' ', OPTIONS.Target(is).Label];
            end
            Comment = [Comment, ', ', OPTIONS.ScoutFunc];
            if ~strcmpi(OPTIONS.ScoutFunc, 'All')
                 Comment = [Comment, ' ' OPTIONS.ScoutTime];
            end
        % Incomplete lists of sources  (not all the sources are present)
        elseif strcmpi(sInput.DataType, 'results') && (length(sInput.RowNames) * sInput.nComponents < nSignals)
            Comment = [Comment, ': ', num2str(length(sInput.RowNames)), ' sources'];
            % Switch the datatype to "scout"
            sInput.DataType = 'scout';
            % Convert source indices to strings
            if ~iscell(sInput.RowNames)
                sInput.RowNames = cellfun(@num2str, num2cell(sInput.RowNames), 'UniformOutput', 0);
            end
        % Single input
        elseif (length(sInput.RowNames) == 1)
            if iscell(sInput.RowNames)
                Comment = [Comment, ': ' sInput.RowNames{1}];
            else
                Comment = [Comment, ': #', num2str(sInput.RowNames(1))];
            end
        end
%         RowNames = sInput.RowNames;
%         sInput.RowNames = cell(length(RowNames)*length(RowNames),1);
%         for i=1:length(RowNames)
%             for j =1:length(RowNames)
%                 sInput.RowNames{j+(i-1)*length(RowNames),1}=[RowNames{i} '-' RowNames{j}];
%             end
%         end
        % ===== SAVE FILE / COMPUTE AVERAGE =====
        % Save each as an independent file
        if ~OPTIONS.isAvgOutput
            nAvg = 1;
            OutputFiles{end+1} = SaveFile(DirectPAC, sInput.iStudy, sInputA(iFile).FileName, sInputRef, sInput, Comment, nAvg, OPTIONS, []);
            %OutputFiles{end+1} = SaveFile(DirectPAC, LowFreqs, HighFreqs, nAvg, sInput.iStudy, sInputA(iFile).FileName, sInput, Comment, OPTIONS);
        else
            % Compute online average of the connectivity matrices
            if isempty(DirectPAC_avg)
                DirectPAC_avg = DirectPAC ./ length(sInputA);
            else
                DirectPAC_avg = DirectPAC_avg + DirectPAC ./ length(sInputA);
            end
            nAvg = nAvg + 1;
        end
        
    end
    
    % ===== SAVE AVERAGE =====
    if OPTIONS.isAvgOutput
        % Output study, in case of average
        [tmp, iOutputStudy] = bst_process('GetOutputStudy', sProcess, sInputA);
        % Save file
        OutputFiles{end+1} = SaveFile(DirectPAC_avg, iOutputStudy, [], sInputRef, sInput, Comment, nAvg, OPTIONS, [])
%         OutputFiles{1} = SaveFile(DirectPAC_avg, LowFreqs, HighFreqs, nAvg, iOutputStudy, [], sInput, Comment, OPTIONS);
    end
end
function TF = ana_vec(F,BandBounds,Fs,width,SegmentOverlap) 
    lower_bin = BandBounds(1);
    upper_bin = BandBounds(2);
    N = size(F,1);
    TF = [];
    w = floor(width*Fs*(1/lower_bin));
    ow = floor(w*SegmentOverlap);
    nF = length(lower_bin:upper_bin);
    for i = 1:N              
        [S,fvec,tvec]=spectrogram([F(i,end:-1:1) F(i,:) F(i,end:-1:1)],w,ow,lower_bin:upper_bin,Fs);
        ind = (tvec > length(F)/Fs) & (tvec <= 2*length(F)/Fs);
        S = S(:,ind);

        if isempty(TF)
            TF = zeros(sum(ind),length(fvec),N);
        end
        TF(:,:,i) = S';
    end
end
function TF = amp_vec(F,BandBounds,Fs,width,TFmethod) 
    lower_bin = BandBounds(1);
    upper_bin = BandBounds(2);
    N = size(F,1);
    T = size(F,2);
    TF = [];
    cf = (lower_bin + floor((upper_bin- lower_bin)/2));
    w = floor(width*Fs*(1/lower_bin));
    for i = 1:N
        if strcmp(TFmethod,'wavelet')
            F1 = ampvec(cf, F(i,:)', Fs, width)';          
        elseif strcmp(TFmethod,'hilbert')  
            % Band-pass filter in one frequency band
            Fband = process_bandpass('Compute', F(i,:), Fs, lower_bin, upper_bin, [], 1);
            %Fband = Fband(T+(1:T));
            % Apply Hilbert transform
            F1 = abs(hilbert(Fband')');
            F1 = F1(ceil(0.1*T)+(1:floor(0.8*T))-1);
        elseif strcmp(TFmethod,'stft')            
            [S,fvec,tvec]=spectrogram([F(i,end:-1:1) F(i,:) F(i,end:-1:1)],w,w-1,lower_bin:upper_bin,Fs);
            ind = (tvec > length(F)/Fs) & (tvec <= 2*length(F)/Fs);
            F1 = mean(abs(S(:,ind)));
        end
        if isempty(TF)
            TF = zeros(N,length(F1));
        end
        TF(i,:) = F1;
    end
end
function TF = bp_vec(F,BandBounds,Fs,width,TFmethod)
    N = size(F,1);
    T = size(F,2);
    TF = [];
    lower_bin = BandBounds(1);
    upper_bin = BandBounds(2);
        
    for i = 1:N
        if strcmp(TFmethod,'wavelet')
            F1 = bpvec((lower_bin + floor((upper_bin- lower_bin)/2)),F(i,:)', Fs, width);

        elseif strcmp(TFmethod,'hilbert') || strcmp(TFmethod,'stft') 
            F1 = process_bandpass('Compute', F(i,:), Fs, lower_bin, upper_bin, [], 1)';
            if strcmp(TFmethod,'hilbert')
                F1 = F1(ceil(0.1*T)+(1:floor(0.8*T))-1);
            end
        end
        if isempty(TF)
            TF = zeros(N,length(F1));
        end
        TF(i,:) = F1';

    end
    
end
function TF = ph_vec(F,BandBounds,Fs,width,TFmethod) 

    lower_bin = BandBounds(1);
    upper_bin = BandBounds(2);
    N = size(F,1);
    T = size(F,2);
    TF = [];
    cf = (lower_bin + floor((upper_bin- lower_bin)/2));
    w = floor(width*Fs*(1/lower_bin));
    for i = 1:N
        if strcmp(TFmethod,'wavelet')
            F1 = phasevec(cf, F(i,:)', Fs, width)';          
        elseif strcmp(TFmethod,'hilbert')  
            % Band-pass filter in one frequency band
            Fband = process_bandpass('Compute', F(i,:), Fs, lower_bin, upper_bin, [], 1);
            % Apply Hilbert transform
            F1 = angle(hilbert(Fband')');
            F1 = F1(ceil(0.1*T)+(1:floor(0.8*T))-1);
        elseif strcmp(TFmethod,'stft')            
            [S,fvec,tvec]=spectrogram([F(i,end:-1:1) F(i,:) F(i,end:-1:1)],w,w-1,lower_bin:upper_bin,Fs);
            ind = (tvec > length(F)/Fs) & (tvec <= 2*length(F)/Fs);
            F1 = mean(angle(S(:,ind)));
        end
        if isempty(TF)
            TF = zeros(N,length(F1));
        end
        TF(i,:) = F1;
    end


end
function mival = mi_measure(phase_sig, amp_sig)
% function mival = mi_measure(phase_sig, amp_sig)
%
% Returns a value for the MI measure calculated between two signals.
% (Functionality to deal with multiple trials will be added soon)
%
% INPUTS:
%
% phase_sig - the instantaneous phase values for a signal which has been
% filtered for a lower, modulating frequency, passed as a column vector
%
% amp_sig - the amplitude values for a signal which has been filtered for a
% higher, modulated frequency, passed as a column vector 
%
% Author: Angela Onslow, May 2010

    num_compX = size(phase_sig, 1);
    num_compY = size(amp_sig, 1);
    mival = zeros(num_compX,num_compY);
    for countX = 1:num_compX
        for countY = 1:num_compY
            %Create composite signal
            z = amp_sig(countY,:).*exp(1i*phase_sig(countX,:));
            m_raw= mean(z);  %Compute the mean length of composite signal.
            mival(countX,countY) =  abs((m_raw));
        end
    end

end

function coh = coherence(sigX, sigY, bandBounds, sRate, width, SegmentOverlap)
    
    TFX = ana_vec(sigX,bandBounds,sRate,width,SegmentOverlap); 
    TFY = ana_vec(sigY,bandBounds,sRate,width,SegmentOverlap); 
    coh = 0;
    for f = 1:size(TFX,2)
        scoh = abs(sum(TFX(:,f,1).*conj(TFY(:,f,1))))^2/norm(TFX(:,f,1))^2/norm(TFY(:,f,1))^2;
        coh = coh + scoh;
    end
    coh = coh / size(TFX,2);
    
%     [Gxy, pValues, freq] = bst_cohn(sigX, sigY, sRate, [], SegmentOverlap, [], 0, [], 0);
%     coh = zeros(num_compX,num_compY);
%     nf = 0;
%     for f = 1:length(freq)
%         if freq(f) > bandBounds(1) && freq(f) < bandBounds(2)
%             coh = coh + Gxy(:,:,f);
%             nf = nf + 1;
%         end
%     end
%     coh = coh / nf;
end
%% ===== SAVE FILE =====
function NewFile = SaveFile(R, iOuptutStudy, DataFile, sInputA, sInputB, Comment, nAvg, OPTIONS, FreqBands)
    NewFile = [];
    bst_progress('text', 'Saving results...');

    % ===== PREPARE OUTPUT STRUCTURE =====
    % Create file structure
    FileMat = db_template('timefreqmat');
    FileMat.TF        = R;
    FileMat.Comment   = Comment;
    FileMat.DataType  = sInputB(1).DataType;
    FileMat.Freqs     = OPTIONS.Freqs;
    FileMat.Method    = OPTIONS.Method;
    FileMat.DataFile  = file_win2unix(DataFile);
    FileMat.nAvg      = nAvg;
    % Time vector
    if strcmpi(OPTIONS.Method, 'plvt')
        FileMat.Time      = sInputB.Time;
        FileMat.TimeBands = [];
    else
        FileMat.Time      = sInputB.Time([1,end]);
        FileMat.TimeBands = {OPTIONS.Method, sInputB.Time(1), sInputB.Time(end)};
    end
    % Measure
    if strcmpi(OPTIONS.Method, 'plv') || strcmpi(OPTIONS.Method, 'plvt')
        FileMat.Measure   = 'none';
    else
        FileMat.Measure   = 'other';
    end
    % Row names: NxM
    FileMat.RefRowNames = sInputA.RowNames;
    FileMat.RowNames    = sInputB.RowNames;
    % Atlas 
    if isstruct(OPTIONS.Target)
        % Save the atlas in the file
        FileMat.Atlas = db_template('atlas');
        FileMat.Atlas.Name   = OPTIONS.ProcessName;
        FileMat.Atlas.Scouts = OPTIONS.Target;
    elseif ~isempty(sInputB.Atlas)
        FileMat.Atlas = sInputB.Atlas;
    end
    if ~isempty(sInputB.SurfaceFile)
        FileMat.SurfaceFile = sInputB.SurfaceFile;
    end
    if ~isempty(sInputB.GridLoc)
        FileMat.GridLoc = sInputB.GridLoc;
    end
    % History: Computation
    FileMat = bst_history('add', FileMat, 'compute', ['Connectivity measure: ', OPTIONS.Method, ' (see the field "Options" for input parameters)']);
    % Save options structure
    FileMat.Options = OPTIONS;
    % Apply time and frequency bands
    if ~isempty(FreqBands)
        FileMat = process_tf_bands('Compute', FileMat, FreqBands, []);
        if isempty(FileMat)
            bst_report('Error', OPTIONS.ProcessName, [], 'Error computing the frequency bands.');
            return;
        end
    end
    
    % ===== PROCESS SCOUTS =====
    % Process scouts: call aggregating function
    if (OPTIONS.isScout) && strcmpi(OPTIONS.ScoutTime, 'after') && ~strcmpi(OPTIONS.ScoutFunc, 'all')
        
        sScoutsA = OPTIONS.Target;

        sScoutsB = [];
        
        FileMat = process_average_rows('ProcessConnectScouts', FileMat, OPTIONS.ScoutFunc, sScoutsA, sScoutsB);
    end
    
    % ===== OPTIMIZE STORAGE FOR SYMMETRIC MATRIX =====
    % Keep only the values below the diagonal
    if FileMat.Options.isSymmetric && (size(FileMat.TF,1) == length(FileMat.RowNames)^2)
        FileMat.TF = process_compress_sym('Compress', FileMat.TF);
    end
        
    % ===== SAVE FILE =====
    % Get output study
    sOutputStudy = bst_get('Study', iOuptutStudy);
    % File tag
    if (length(FileMat.RefRowNames) == 1)
        fileTag = 'connect1';
    else
        fileTag = 'connectn';
    end
    % Output filename
    NewFile = bst_process('GetNewFilename', bst_fileparts(sOutputStudy.FileName), ['timefreq_' fileTag '_' OPTIONS.Method]);
    % Save file
    bst_save(NewFile, FileMat, 'v6');
    % Add file to database structure
    db_add_data(iOuptutStudy, NewFile, FileMat);
end
function [c, lags] = lagged_corr(X,Y,maxlag)
    
    N = length(X);
    if length(Y) ~=N
        c=0; lags=0;
        return;
    end
    % if val < 3
    %     maxlag = N;
    % end


    if maxlag == 0
        z = corrcoef(X,Y);
        c = z(1,2);
        lags = 0;
        return;
    end
    K = N - 2*maxlag;
    ms = -maxlag:10:maxlag;
    for t=1:length(ms)
        m=ms(t);
        z = corrcoef(X(m+maxlag+(1:K)),Y(maxlag+(1:K)-1));
        c(t) = z(1,2);
    end
    [val,ind]=max(c);
    c=val;
    lags = ms(ind);
end
