function varargout = process_CFC( varargin )
% PROCESS_CFC: Compute the Cross-Frequency Coupling between or within time
% series.
%
% DOCUMENTATION: For more detail, please see the tutorial in https://bsp.hackpad.com/Cross-Frequency-Coupling-cChe95lhDHz

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
% Authors: Hui-Ling Chan 2015

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Cross-Frequency Coupling';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Frequency';
    sProcess.Index       = 660;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw',      'data',     'results',  'matrix'};
    sProcess.OutputTypes = {'timefreq', 'timefreq', 'timefreq', 'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % ==== INPUT ====
    sProcess.options.label_in.Comment = '<HTML><B><U>Input options</U></B>:';
    sProcess.options.label_in.Type    = 'label';
    sProcess.options.result_comm.Comment    = 'Comment: ';
    sProcess.options.result_comm.Type       = 'text';
    sProcess.options.result_comm.Value      = '';
    % === TIME WINDOW
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % === SENSOR SELECTION
    sProcess.options.target_data.Comment    = 'Sensor types or names (empty=all): ';
    sProcess.options.target_data.Type       = 'text';
    sProcess.options.target_data.Value      = 'MEG, EEG';
    sProcess.options.target_data.InputTypes = {'data', 'raw'};
    % === SCOUTS SELECTION
    sProcess.options.scouts.Comment    = 'Use scouts';
    sProcess.options.scouts.Type       = 'scout_confirm';
    sProcess.options.scouts.Value      = {};
    sProcess.options.scouts.InputTypes = {'results'};
    % === SCOUT FUNCTION ===
    sProcess.options.scoutfunc.Comment    = {'Mean', 'Max', 'PCA', 'Std', 'All', 'Scout function:'};
    sProcess.options.scoutfunc.Type       = 'radio_line';
    sProcess.options.scoutfunc.Value      = 1;
    sProcess.options.scoutfunc.InputTypes = {'results'};
    % === SCOUT TIME ===
    sProcess.options.scouttime.Comment    = {'Before', 'After', 'When to apply the scout function:'};
    sProcess.options.scouttime.Type       = 'radio_line';
    sProcess.options.scouttime.Value      = 1;
    sProcess.options.scouttime.InputTypes = {'results'};
    % === ROW NAMES
    sProcess.options.target_tf.Comment    = 'Row names or indices (empty=all): ';
    sProcess.options.target_tf.Type       = 'text';
    sProcess.options.target_tf.Value      = '';
    sProcess.options.target_tf.InputTypes = {'timefreq', 'matrix'};

    % ==== ESTIMATOR ====
    sProcess.options.label_pac.Comment = '<HTML><BR><B><U>Estimator options</U></B>:';
    sProcess.options.label_pac.Type    = 'label';
    % === PAC MEASURES ===
    sProcess.options.pacmeasure.Comment    = {'AEC', 'ESC', 'ESP', 'MI', 'CFC measure:'};
    sProcess.options.pacmeasure.Type       = 'radio_line';
    sProcess.options.pacmeasure.Value      = 1;
    % === TIME LAGGED
    sProcess.options.tlag.Comment = 'Time-lagged';
    sProcess.options.tlag.Type    = 'checkbox';
    sProcess.options.tlag.Value   = 1;
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
    sProcess.options.nestedwidth.Comment = 'Nested frequency step (high):';
    sProcess.options.nestedwidth.Type    = 'value';
    sProcess.options.nestedwidth.Value   = {0.75, 'Hz', 2};   
%     % === TF METHOD  ===
%     sProcess.options.tfmethod.Comment    = {'Hilbert', 'Wavelet', 'STFT', 'TF method:'};
%     sProcess.options.tfmethod.Type       = 'radio_line';
%     sProcess.options.tfmethod.Value      = 1;
  
    %  width - width of the wavelet filter; default = 7
%  nfft - the number of points in fft; default = 200
%  num_shf - the number of shuffled data sets to use during significance
%  testing; default = 50
%  alpha - significance value to use; default = 0.05 

%     % ==== ESTIMATOR ====
%     sProcess.options.label_ad.Comment = '<HTML><BR><B><U>Advanced options</U></B>:';
%     sProcess.options.label_ad.Type    = 'label';
%     % ===  width of the wavelet filter 
%     sProcess.options.width.Comment = 'Width of the wavelet filter (cycles):';
%     sProcess.options.width.Type    = 'value';
%     sProcess.options.width.Value   = {7, '', 0};
    % === the number of points in fft 
%     sProcess.options.nfft.Comment = 'Number of points in fft:';%only for coherence computation
%     sProcess.options.nfft.Type    = 'value';
%     sProcess.options.nfft.Value   = {200, '', 0};   
%     % === the number of shuffled data sets to use during significance 
%     sProcess.options.num_shf.Comment = 'Number of shuffled data sets:';
%     sProcess.options.num_shf.Type    = 'value';
%     sProcess.options.num_shf.Value   = {50, '', 0};
%     % === alpha 
%     sProcess.options.alpha.Comment = 'Significance value to use:';
%     sProcess.options.alpha.Type    = 'value';
%     sProcess.options.alpha.Value   = {0.05, '', 3};
    % === Width 
%     sProcess.options.width.Comment = 'Cycles of wavelet:';
%     sProcess.options.width.Type    = 'value';
%     sProcess.options.width.Value   = {7, 'Cycles', 0}; 
%     % === WINDOW LENGTH
%     sProcess.options.winlength.Comment = 'Estimator window length: ';
%     sProcess.options.winlength.Type    = 'value';
%     sProcess.options.winlength.Value   = {0.128, 'ms ', 1};
%     % === Low bound
%     sProcess.options.winoverlap.Comment = 'Overlap percentage: ';
%     sProcess.options.winoverlap.Type    = 'value';
%     sProcess.options.winoverlap.Value   = {0.75, '% ', 1};

    % ==== OUTPUT ====
    sProcess.options.label_out.Comment = '<HTML><BR><U><B>Output configuration</B></U>:';
    sProcess.options.label_out.Type    = 'label';
    % === SCOUT FUNCTION ===
    sProcess.options.pacmatrix.Comment    = {'Auto', 'Cross', 'Both', 'Output:'};
    sProcess.options.pacmatrix.Type       = 'radio_line';
    sProcess.options.pacmatrix.Value      = 1;
    % === AVERAGE OUTPUT FILES
    sProcess.options.avgoutput.Comment = 'Save average CFC across trials (one output file only)';
    sProcess.options.avgoutput.Type    = 'checkbox';
    sProcess.options.avgoutput.Value   = 1;
    % === SAVE PAC MAPS
    sProcess.options.savemax.Comment = 'Save only the maximum CFC values';
    sProcess.options.savemax.Type    = 'checkbox';
    sProcess.options.savemax.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputA) %#ok<DEFNU>
    OutputFiles = {};
    
    % ===== GET OPTIONS =====
    if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
        OPTIONS.TimeWindow = sProcess.options.timewindow.Value{1};
    else
        OPTIONS.TimeWindow = [];
    end
    result_comment = sProcess.options.result_comm.Value;
    if ~isempty(result_comment)
        result_comment = [result_comment ': '];
    end
    % Get and check frequencies
    OPTIONS.BandNesting = sProcess.options.nesting.Value{1};
    OPTIONS.BandNested  = sProcess.options.nested.Value{1};
    if (min(OPTIONS.BandNesting) < 0.5)
        bst_report('Error', sProcess, [], 'This function cannot be used to estimate PAC for nesting frequencies below 1Hz.');
        return;
    end
%     if (max(OPTIONS.BandNesting) > min(OPTIONS.BandNested))
%         bst_report('Error', sProcess, [], 'The low and high frequency band cannot overlap.');
%         return;
%     end
    % Get target
    if ismember(sInputA(1).FileType, {'data','raw'}) && isfield(sProcess.options, 'target_data') && ~isempty(sProcess.options.target_data.Value)
        OPTIONS.Target = sProcess.options.target_data.Value;
%     elseif strcmpi(sInputA(1).FileType, 'results') && isfield(sProcess.options, 'scouts') && ~isempty(sProcess.options.scouts.Value)
%         OPTIONS.Target = sProcess.options.scouts.Value{2};
    elseif ismember(sInputA(1).FileType, {'timefreq', 'matrix'}) && isfield(sProcess.options, 'target_tf') && ~isempty(sProcess.options.target_tf.Value)
        OPTIONS.Target = sProcess.options.target_tf.Value;
    else
        OPTIONS.Target = [];
    end
    % All other options
    %OPTIONS.MaxSignals   = sProcess.options.max_block_size.Value{1};
    %OPTIONS.isParallel   = sProcess.options.parallel.Value && (exist('matlabpool', 'file') ~= 0);
    %OPTIONS.isMex        = sProcess.options.ismex.Value;
    OPTIONS.isSaveMax    = sProcess.options.savemax.Value;
    OPTIONS.isTimeLag    = sProcess.options.tlag.Value;
    OPTIONS.isAvgOutput  = sProcess.options.avgoutput.Value;
    if (length(sInputA) == 1)
        OPTIONS.isAvgOutput = 0;
    end
    
    switch (sProcess.options.pacmatrix.Value)
        case 1, OPTIONS.PACmatrix = 'auto';
        case 2, OPTIONS.PACmatrix = 'cross';
        case 3, OPTIONS.PACmatrix = 'both';
    end
    
    switch (sProcess.options.pacmeasure.Value)
        case 1, OPTIONS.PACmeasure = 'aec';
        case 2, OPTIONS.PACmeasure = 'esc';
        case 3, OPTIONS.PACmeasure = 'esp';
        case 4, OPTIONS.PACmeasure = 'mi';
    end
    sProcess.options.tfmethod.Value = 1;
    switch (sProcess.options.tfmethod.Value)
        case 1, OPTIONS.TFmethod = 'hilbert';
        case 2, OPTIONS.TFmethod = 'wavelet';
        case 3, OPTIONS.TFmethod = 'stft';
    end
    OPTIONS.nfft = 256;
    OPTIONS.Width = 7;
%     OPTIONS.width = sProcess.options.width.Value{1};
%     OPTIONS.nfft = sProcess.options.nfft.Value{1};
%     OPTIONS.WinOverlap = sProcess.options.winoverlap.Value{1};
%     OPTIONS.SegmentLength = sProcess.options.winlength.Value{1};
%    OPTIONS.Width = sProcess.options.width.Value{1};
%     OPTIONS.num_shf = sProcess.options.num_shf.Value{1};
%     OPTIONS.alpha = sProcess.options.alpha.Value{1};
    % ===== GET SCOUTS OPTIONS =====
    if strcmpi(sInputA(1).FileType, 'results') && isfield(sProcess.options, 'scouts') && isfield(sProcess.options.scouts, 'Value')
        % Selected scouts
        sScouts = sProcess.options.scouts.Value{2};
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
        %if strcmpi(OPTIONS.ScoutTime, 'before') && ismember(OPTIONS.ScoutFunc, {'max', 'std'})
        %    bst_report('Error', sProcess, [], 'Scout functions MAX and STD should not be applied before estimating the PAC.');
        %    return;
        %end
        if strcmpi(OPTIONS.ScoutTime, 'after') && strcmpi(OPTIONS.ScoutFunc, 'pca')
            bst_report('Error', sProcess, [], 'Scout function PCA cannot be applied after estimating the PAC.');
            return;
        end
%         % Set input/output scouts functions
%         if ~isempty(sScouts)
%             OPTIONS.Target = sScouts;
%             % Apply function before: get all the scouts time series in advance
%             if strcmpi(OPTIONS.ScoutTime, 'before')
%                 [OPTIONS.Target.Function] = deal(OPTIONS.ScoutFunc);
%             % Apply function after: Get all the time series of all the scouts
%             elseif strcmpi(OPTIONS.ScoutTime, 'after')
%                 [OPTIONS.Target.Function] = deal('All');
%             end
%         end
         % Selected scouts
        AtlasList = sProcess.options.scouts.Value;
        % Set input/output scouts functions
        if ~isempty(AtlasList)
            OPTIONS.Target = AtlasList;
            % Apply function before: get all the scouts time series in advance
            if strcmpi(OPTIONS.ScoutTime, 'before')
                LoadOptions.TargetFunc = OPTIONS.ScoutFunc;
            % Apply function after: Get all the time series of all the scouts
            elseif strcmpi(OPTIONS.ScoutTime, 'after')
                LoadOptions.TargetFunc = 'all';
            end
        end
    end
    
    % ===== INITIALIZE =====
    % Initialize output variables
    DirectPAC_avg = [];
    LowFreqs  = [];
    HighFreqs = [];
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
    
    BandNesting_vec = OPTIONS.BandNesting(1):sProcess.options.nestingwidth.Value{1}:OPTIONS.BandNesting(2);
    BandNested_vec  = OPTIONS.BandNested(1):sProcess.options.nestedwidth.Value{1}:OPTIONS.BandNested(2);
    
    % Loop over input files
    for iFile = 1:length(sInputA)
        % ===== LOAD SIGNALS =====
        bst_progress('text', sprintf('PAC: Loading input file (%d/%d)...', iFile, length(sInputA)));
        bst_progress('set', round(startValue + (iFile-1) / length(sInputA) * 100));
        % Load input signals 
        [sInput, nSignals, iRows] = bst_process('LoadInputFile', sInputA(iFile).FileName, OPTIONS.Target, OPTIONS.TimeWindow, LoadOptions);
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

        
        % Definitions
%         bandNesting = OPTIONS.BandNesting;
%         bandNested = OPTIONS.BandNested;
%         fmin = min(bandNesting);
%         fmax = sRate/3;
%         numfreqs = round(sRate/9);
%         fstep = sProcess.options.nestingwidth.Value{1};
%         % Calculate center frequencies
%         temp1 = (0:numfreqs-1) * fstep;
%         temp2 = logspace(log10(fmin), log10(fmax), numfreqs);
%         temp2 = (temp2-temp2(1)) * ((temp2(end)-temp1(end)) / temp2(end)) + temp2(1);
%         chirpCenterFreqs = temp1 + temp2;
%         % Remove unused frequencies
%         chirpCenterFreqs(chirpCenterFreqs > max(bandNested)) = [];      %%% ESTHER
%         chirpCenterFreqs((chirpCenterFreqs < min(bandNested)) & (chirpCenterFreqs >= max(bandNesting))) = [];      %%% ESTHER
% % %         % Indices of center frequencies in the upper frequency range
% % %         hfreq = find( chirpCenterFreqs >= min(bandNested) );
% % %         % Number of cf bins to evaluate for PAC with lower-frequency oscillations
% % %         % lfreq = find(chirpCenterFreqs < min(bandNested));
% % %         lfreq = find(chirpCenterFreqs < max(bandNesting));   %%% ESTHER
%         BandNested_vec = chirpCenterFreqs(chirpCenterFreqs >= min(bandNested) );
%         BandNesting_vec = chirpCenterFreqs(chirpCenterFreqs < max(bandNesting) );
    
        % ===== COMPUTE PAC MEASURE =====
        % Number of blocks of signals
        %MAX_BLOCK_SIZE = 1;%OPTIONS.MaxSignals;
        %nBlocks = ceil(nSignals / MAX_BLOCK_SIZE);
        DirectPAC = [];
        % Display processing time
        %disp(sprintf('BST> PAC: Processing %d blocks of %d signals each.', nBlocks, MAX_BLOCK_SIZE));
        % Process each block of signals
        %for iBlock = 1:nBlocks
            tic
%             bst_progress('text', sprintf('PAC: File %d/%d - Block %d/%d', iFile, length(sInputA), iBlock, nBlocks));
%             bst_progress('set', round(startValue + (iFile-1)/length(sInputA)*100 + iBlock/nBlocks*100));    
%             % Indices of the signals
%             iSignals = (iBlock-1)*MAX_BLOCK_SIZE+1 : min(iBlock*MAX_BLOCK_SIZE, nSignals);
            % Get signals to process
            if ~isempty(sInput.ImagingKernel)
                Fblock = sInput.ImagingKernel * sInput.Data;
            else
                Fblock = sInput.Data;
            end
            
            % ===== APPLY SOURCE ORIENTATION =====
            if strcmpi(sInput.DataType, 'results')
                % Number of values per vertex
                switch (sInput.nComponents)
                    case 0
                        error('PAC metrics are not supported for mixed source models.');
                    case 1
                        % Nothing to do
                    case 2
                        Fblock = (Fblock(1:2:end,:) + Fblock(2:2:end,:)) / 2;
                        sInput.RowNames = sInput.RowNames(1:2:end);
                    case 3
                        Fblock = (Fblock(1:3:end,:) + Fblock(2:3:end,:) + Fblock(3:3:end,:)) / 3;
                        sInput.RowNames = sInput.RowNames(1:3:end);
                end
            end
            %[DirectPAC_block, LowFreqs, HighFreqs] = bst_pac(Fblock, sRate, OPTIONS.BandNesting, OPTIONS.BandNested, OPTIONS.isParallel, OPTIONS.isMex);
            
            % [pacmat, freqvec_ph, freqvec_amp] = find_pac_shf (sig_pac, Fs, measure, ...
            %     sig_mod, ph_freq_vec, amp_freq_vec, plt, waitbar, width, nfft, num_shf, alpha,...
            %     dataname, sig_pac_name, sig_mod_name)

            n=0;
            for iSigX = 1:nSignals
                for iSigY = 1:nSignals
                    if iSigX == iSigY && strcmpi(OPTIONS.PACmatrix, 'cross')
                        continue;
                    end
                    if iSigX ~= iSigY && strcmpi(OPTIONS.PACmatrix, 'auto')
                        continue;
                    end                           
%                     OPTIONS.num_shf = 0;
%                     OPTIONS.alpha = 1;
%                     [pacmat, LowFreqs, HighFreqs] = find_pac_shf(Fblock(iSigY,:), sRate, OPTIONS.PACmeasure, Fblock(iSigX,:), BandNesting_vec,BandNested_vec,'n',1,OPTIONS.Width, OPTIONS.nfft, OPTIONS.num_shf, OPTIONS.alpha);
                    [pacmat, LowFreqs, HighFreqs] = find_pac_var_bandwidth(Fblock(iSigY,:), sRate, OPTIONS.PACmeasure, Fblock(iSigX,:), BandNesting_vec, BandNested_vec, OPTIONS.Width, OPTIONS.nfft, OPTIONS.TFmethod, OPTIONS.isTimeLag);
                    if isempty(DirectPAC)
                        RowNames = sInput.RowNames;
                        if  strcmpi(OPTIONS.PACmatrix, 'cross')
                            sInput.RowNames = cell((nSignals^2)-nSignals,1);
                            DirectPAC =  zeros((nSignals^2)-nSignals,1,size(pacmat,2),size(pacmat,1));
                        elseif strcmpi(OPTIONS.PACmatrix, 'auto') 
                            DirectPAC =  zeros(nSignals,1,size(pacmat,2),size(pacmat,1));
                        elseif strcmpi(OPTIONS.PACmatrix, 'both')
                            sInput.RowNames = cell(length(RowNames)^2,1);
                            DirectPAC =  zeros(nSignals*nSignals,1,size(pacmat,2),size(pacmat,1));
                        end
                    end
                    n=n+1;
                    DirectPAC(n,1,:,:) = pacmat'; 
                    if ~strcmpi(OPTIONS.PACmatrix, 'auto')
                        sInput.RowNames{n,1}=[RowNames{iSigX} '-' RowNames{iSigY}];
                    end

                end
            end

            LowFreqs = LowFreqs';
            HighFreqs = HighFreqs';
            
        % ===== APPLY SOURCE ORIENTATION =====
        % Unconstrained sources => SUM for each point
        if ismember(sInput.DataType, {'results','scout','matrix'}) && ~isempty(sInput.nComponents) && (sInput.nComponents ~= 1)
            [DirectPAC, sInput.GridAtlas, sInput.RowNames] = bst_source_orient([], sInput.nComponents, sInput.GridAtlas, DirectPAC, 'mean', sInput.DataType, sInput.RowNames);
        end
        % ===== PROCESS SCOUTS =====
        % Get scouts
        isScout = ~isempty(OPTIONS.Target) && (isstruct(OPTIONS.Target) || iscell(OPTIONS.Target)) && isfield(sInput, 'Atlas') && isfield(sInput.Atlas, 'Scouts') && ~isempty(sInput.Atlas.Scouts);    
        if isScout
            sScouts = sInput.Atlas.Scouts;
        end
        % If the scout function has to be applied AFTER the PAC computation
        if isScout && strcmpi(OPTIONS.ScoutTime, 'after') && ~strcmpi(OPTIONS.ScoutFunc, 'all')
            nScouts = length(sScouts);
            DirectPAC_scouts = zeros(nScouts, size(DirectPAC,2), size(DirectPAC,3), size(DirectPAC,4));
            iVerticesAll = [1, cumsum(cellfun(@length, {sScouts.Vertices})) + 1];
            % For each unique row name: compute a measure over the clusters values
            for iScout = 1:nScouts
                iScoutVert = iVerticesAll(iScout):iVerticesAll(iScout+1)-1;
                F = reshape(DirectPAC(iScoutVert,:,:,:), length(iScoutVert), []);
                F = bst_scout_value(F, OPTIONS.ScoutFunc);
                DirectPAC_scouts(iScout,:,:,:) = reshape(F, [1, size(DirectPAC,2), size(DirectPAC,3), size(DirectPAC,4)]);
            end
            % Save only the requested rows
            sInput.RowNames = {sScouts.Label};
            DirectPAC = DirectPAC_scouts;
        end
        % ===== FILE COMMENT =====
        % Base comment
        Comment = result_comment;
        if OPTIONS.isTimeLag
            Comment = [Comment 'Lagged'];
        end
        if OPTIONS.isSaveMax
            Comment = [Comment 'Max' upper(OPTIONS.PACmeasure) '(' upper(OPTIONS.TFmethod) ')'];
        else
            Comment = [Comment upper(OPTIONS.PACmeasure) '(' upper(OPTIONS.TFmethod) ')'];           
        end

        % Time window (RAW only)
        if ~isempty(strfind(sInputA(iFile).Condition, '@raw'))
            Comment = [Comment, sprintf('(%ds-%ds)', round(OPTIONS.TimeWindow))];
        end
        % Scouts
        if isScout && (length(sScouts) < 6)
            Comment = [Comment, ':'];
            for is = 1:length(sScouts)
                Comment = [Comment, ' ', sScouts(is).Label];
            end
            Comment = [Comment, ', ', OPTIONS.ScoutFunc];
            if ~strcmpi(OPTIONS.ScoutFunc, 'All')
                 Comment = [Comment, ' ' OPTIONS.ScoutTime];
            end
        % Single input
        elseif (length(sInput.RowNames) == 1)
            if iscell(sInput.RowNames)
                Comment = [Comment, ': ' sInput.RowNames{1}];
            else
                Comment = [Comment, ': #', num2str(sInput.RowNames(1))];
            end
        end
        
        % ===== SAVE FILE / COMPUTE AVERAGE =====
        % Save each as an independent file
        if ~OPTIONS.isAvgOutput
            nAvg = 1;
            OutputFiles{end+1} = SaveFile(DirectPAC, LowFreqs, HighFreqs, nAvg, sInput.iStudy, sInputA(iFile).FileName, sInput, Comment, OPTIONS);
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
        OutputFiles{1} = SaveFile(DirectPAC_avg, LowFreqs, HighFreqs, nAvg, iOutputStudy, [], sInput, Comment, OPTIONS);
    end
end


%% ========================================================================
%  ===== SUPPORT FUNCTIONS ================================================
%  ========================================================================

%% ===== SAVE FILE =====
function NewFile = SaveFile(DirectPAC, LowFreqs, HighFreqs, nAvg, iOuptutStudy, DataFile, sInput, Comment, OPTIONS)
    % ===== COMPUTE MAXPAC ======
    % Save directPAC values in returned structure only if requested
    if OPTIONS.isSaveMax
        sPAC.DirectPAC = [];
    else        
        sPAC.DirectPAC = DirectPAC;
    end
    % Get the maximum DirectPAC value for each signal
    [sPAC.ValPAC, indmax] = max(reshape(DirectPAC, size(DirectPAC,1), []), [], 2);
    % Find the pair of low/high frequencies for this maximum
    [imaxl, imaxh]   = ind2sub([size(DirectPAC,3), size(DirectPAC,4)], indmax);
    sPAC.NestingFreq = LowFreqs(imaxl)';
    sPAC.NestedFreq  = HighFreqs(imaxh)';
    % Copy list of frequencies 
    sPAC.LowFreqs  = LowFreqs;
    sPAC.HighFreqs = HighFreqs;

    % ===== PREPARE OUTPUT STRUCTURE =====
    % Create file structure
    FileMat = db_template('timefreqmat');
    FileMat.TF        = sPAC.ValPAC;
    FileMat.Comment   = Comment;
    FileMat.DataType  = sInput.DataType;
    FileMat.RowNames  = sInput.RowNames;
    FileMat.Time      = sInput.Time([1,end]);
    FileMat.Method    = 'pac';
    FileMat.Measure   = 'maxpac';
    FileMat.DataFile  = file_win2unix(DataFile);
    FileMat.nAvg      = nAvg;
    FileMat.Freqs     = 0;
    % Atlas 
    if ~isempty(sInput.Atlas)
        FileMat.Atlas = sInput.Atlas;
    end
    if ~isempty(sInput.GridLoc)
        FileMat.GridLoc = sInput.GridLoc;
    end
    if ~isempty(sInput.GridAtlas)
        FileMat.GridAtlas = sInput.GridAtlas;
    end
    if ~isempty(sInput.SurfaceFile)
        FileMat.SurfaceFile = sInput.SurfaceFile;
    end
    % History: Computation
    FileMat = bst_history('add', FileMat, 'compute', 'PAC measure (see the field "Options" for input parameters)');
    % All the PAC fields and options
    FileMat.Options = OPTIONS;
    FileMat.sPAC    = rmfield(sPAC, 'ValPAC');
    
    % ===== SAVE FILE =====
    % Get output study
    sOutputStudy = bst_get('Study', iOuptutStudy);
    % File tag
    if OPTIONS.isSaveMax
        fileTag = 'timefreq_pac';
    else
        fileTag = 'timefreq_pac_fullmaps';
    end
    % Output filename
    NewFile = bst_process('GetNewFilename', bst_fileparts(sOutputStudy.FileName), fileTag);
    % Save file
    bst_save(NewFile, FileMat, 'v6');
    % Add file to database structure
    db_add_data(iOuptutStudy, NewFile, FileMat);
end

function [pacmat, freqvec_ph, freqvec_amp] = find_pac_var_bandwidth(sig_pac, Fs, measure, ...
    sig_mod, ph_freq_vec, amp_freq_vec, width, nfft, TFmethod, isTLag, num_shf, alpha)

% Checks of input variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 11
    num_shf=0;
end
if nargin < 12
    alpha=0.05;
end
% Check data is columnwise
if size(sig_pac,1)<size(sig_pac,2)
    sig_pac = sig_pac';
end

if size(sig_mod,1)<size(sig_mod,2)
    sig_mod = sig_mod';
end

if (size(sig_pac,2) ~= size(sig_mod,2))
    sprintf('Error - Signals must have the same number of trials')
    return
end

total_num_dp = size(sig_pac,1);
num_trials = 1:size(sig_pac,2);

phfreq_low = min(ph_freq_vec);
phfreq_high = max(ph_freq_vec);
phfreq_bw = diff(ph_freq_vec(1:2));
ampfreq_low = min(amp_freq_vec);
ampfreq_high = max(amp_freq_vec);
ampfreq_bw = diff(amp_freq_vec(1:2));

xbins = ceil((phfreq_high - phfreq_low)/phfreq_bw);
ybins = ceil((ampfreq_high - ampfreq_low)/ampfreq_bw);
alpha = alpha/(xbins*ybins); % Uncomment to use Bonferonni Correction

freqvec_amp = zeros(1,ybins);
freqvec_ph = zeros(1,xbins);
amp_filt_signals = zeros(floor(total_num_dp*0.8),max(num_trials));
ph_filt_signals = zeros(floor(total_num_dp*0.8),max(num_trials));
for y=1:ybins
    freqvec_amp(y) = ampfreq_low+(y-1)*ampfreq_bw;  
end
for x=1:xbins
    freqvec_ph(x) = phfreq_low+(x-1)*phfreq_bw;
end
counter = 0;
maxlag = 0;
countermax = xbins*ybins;
%fprintf('\nCalculating CFC values\n');
if (strcmp(measure, 'esc') || strcmp(measure, 'mi') || strcmp(measure, 'esp') || strcmp(measure, 'aec') )
    pacmat = zeros(ybins,xbins);
    for x=1:xbins

        BandBoundsX(1) = max(0,freqvec_ph(x) - 0.5 * phfreq_bw);
        BandBoundsX(2) = BandBoundsX(1) + phfreq_bw;

        if strcmp(measure, 'esc')
            for i3=1:max(num_trials)
                ph_filt_signals(:,i3) = bp_vec(sig_mod(:,i3)', BandBoundsX, Fs, width, TFmethod)';
            end
        elseif strcmp(measure, 'aec')
            for i3=1:max(num_trials)
                ph_filt_signals(:,i3) = amp_vec(sig_mod(:,i3)', BandBoundsX, Fs, width, TFmethod)';
            end
        else
            for i3=1:max(num_trials)
                ph_filt_signals(:,i3) = ph_vec(sig_mod(:,i3)', BandBoundsX, Fs, width, TFmethod)';
            end
        end        
        
        for s = 1:num_shf
            shuffled_sig_ph{s} = shuffle_esc(ph_filt_signals);
        end
        
        for y=1:ybins
            
            %BandBoundsY(1) = max(0,freqvec_amp(y) - freqvec_ph(x) - phfreq_bw - 1);
            BandBoundsY(1) = max(0,freqvec_amp(y) - 0.5*ampfreq_bw);
            
            if BandBoundsY(1) > freqvec_ph(x)
                %BandBoundsY(2) = freqvec_amp(y) + freqvec_ph(x) + phfreq_bw + 1;
                BandBoundsY(2) = freqvec_amp(y) + ampfreq_bw;
                for i3=1:max(num_trials)
                    amp_filt_signals(:,i3) = amp_vec(sig_pac(:,i3)', BandBoundsY, Fs, width, TFmethod)';
                end
                
                if strcmp(measure, 'esc') || strcmp(measure, 'esp')
                    
                    if isTLag == 1                        
                        maxlag = min(floor(Fs/freqvec_ph(x)/2),floor(total_num_dp*0.8)/4);
                    end
                    %pacmat(y,x) = esc_measure(ph_filt_signals, amp_filt_signals, 'y');
                    pacmat(y,x) = lagged_corr(ph_filt_signals,amp_filt_signals,maxlag);
                    for s = 1:num_shf
                        shf_pacmat_final(s,y,x) = lagged_corr(shuffled_sig_ph{s}, amp_filt_signals, maxlag);
                    end
                else
                    
                    
                    pacmat(y,x) = mi_measure(ph_filt_signals, amp_filt_signals);
                    for s = 1:num_shf
                        shf_pacmat_final(s,y,x) = mi_measure(shuffled_sig_ph{s}, amp_filt_signals);
                    end
                    
                end

            end 
            
            counter = counter+1;
%             if counter == 1
%                 fprintf('%03i%% ', floor((counter/countermax)*100));
%             else
%                 fprintf('\b\b\b\b\b%03i%% ', floor((counter/countermax)*100));
%             end
%             if counter == countermax
%                 fprintf('\n');
%             end
        end
    end
end

%Find mean and standard deviation of shuffled data sets
if num_shf ~= 0
    if strcmp(measure, 'mi')
        for i =1:ybins
            for j=1:xbins
                [shf_data_mean(i,j), shf_data_std(i,j)] = normfit(shf_pacmat_final(:,i,j));
            end
        end

    else
        shf_data_mean = squeeze (mean (shf_pacmat_final, 1));
        shf_data_std = squeeze (std (shf_pacmat_final, 1));
    end
end
% Compute significance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if num_shf ~= 0 
for i = 1:size(pacmat,1)
    for j = 1:size(pacmat,2)
        [h, p] = my_sig_test(pacmat(i,j), squeeze(shf_pacmat_final(:,i,j)), alpha);
        if h == 0
            pacmat(i,j) = 0;
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqvec_amp = ceil(freqvec_amp)';
freqvec_ph = ceil(freqvec_ph)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pacmat, freqvec_ph, freqvec_amp] = find_pac(sig_pac, Fs, measure, ...
    sig_mod, ph_freq_vec, amp_freq_vec, width, nfft, TFmethod, num_shf, alpha)

% Checks of input variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 10
    num_shf=50;
end
if nargin < 11
    alpha=0.05;
end
% Check data is columnwise
if size(sig_pac,1)<size(sig_pac,2)
    sig_pac = sig_pac';
end

if size(sig_mod,1)<size(sig_mod,2)
    sig_mod = sig_mod';
end

if (size(sig_pac,2) ~= size(sig_mod,2))
    sprintf('Error - Signals must have the same number of trials')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up some parameters for clarity
xbins = ceil((max(ph_freq_vec) - min(ph_freq_vec))/(diff(ph_freq_vec(1:2))));
ybins = ceil((max(amp_freq_vec) - min(amp_freq_vec))/(diff(amp_freq_vec(1:2))));
%alpha = alpha/(xbins*ybins); % Uncomment to use Bonferonni Correction

% Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each cell array element of ph_filt_signals and amp_filt_signals has the
% same dimensions as the original signals i.e. number of columns = number of
% trials
if (strcmp(measure, 'esc')) ||(strcmp(measure, 'mi')) || (strcmp(measure, 'esp'))
[filt_sig_mod, filt_sig_pac] = filt_signals(sig_pac, sig_mod, Fs, ...
    ph_freq_vec, amp_freq_vec, measure, width, TFmethod);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create shuffled datasets and distribution of PAC values %%%%%%%%%%%%%%%%%
if num_shf ~= 0
for s = 1:num_shf
    
    if strcmp(measure, 'esc') || strcmp(measure, 'esp')
           
        shuffled_sig_amp = shuffle_esc(filt_sig_pac, Fs);
        shf_pacmat_final(s,:,:) = find_pac_nofilt(shuffled_sig_amp, Fs, 'esc', filt_sig_mod, ph_freq_vec, amp_freq_vec,'n');
        
    end
     

    if strcmp(measure, 'mi')
    
        shuffled_sig_amp = shuffle_esc(filt_sig_pac, Fs);
        shf_pacmat_final(s,:,:) = find_pac_nofilt(shuffled_sig_amp, Fs, measure, filt_sig_mod, ph_freq_vec, amp_freq_vec,'n');
        
    end
    
    if strcmp(measure, 'cfc')
          
        shuffled_sig1 = shuffle_esc(sig_pac, Fs);
        shf_pacmat_final(s,:,:) = find_pac_nofilt(shuffled_sig1, Fs,measure, sig_mod, ph_freq_vec, amp_freq_vec,'n', 0, width, nfft);
        
    end
    
    % Display current computational step to user
    %if waitbar == 1
        if s == 1
            fprintf('%03i%% ', floor((s/num_shf)*100));
        else
            fprintf('\b\b\b\b\b%03i%% ', floor((s/num_shf)*100));
        end
        if s == num_shf
            fprintf('\n');
        end
    %end
end

%Find mean and standard deviation of shuffled data sets
if strcmp(measure, 'mi')
    for i =1:ybins
        for j=1:xbins
            [shf_data_mean(i,j), shf_data_std(i,j)] = normfit(shf_pacmat_final(:,i,j));
        end
    end
    
else
    shf_data_mean = squeeze (mean (shf_pacmat_final, 1));
    shf_data_std = squeeze (std (shf_pacmat_final, 1));
end
end


% Calculate PAC measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(measure, 'esc')  || strcmp(measure, 'esp') 
    [pacmat, freqvec_ph, freqvec_amp] = find_pac_nofilt(filt_sig_pac, Fs, 'esc', filt_sig_mod, ph_freq_vec, amp_freq_vec, 'n', 1);
end

if strcmp(measure, 'mi')
    [pacmat, freqvec_ph, freqvec_amp] = find_pac_nofilt(filt_sig_pac, Fs, measure, filt_sig_mod, ph_freq_vec, amp_freq_vec, 'n', 1);
end

if strcmp(measure, 'cfc')
    [pacmat, freqvec_ph, freqvec_amp] = find_pac_nofilt(sig_pac, Fs, measure, sig_mod, ph_freq_vec, amp_freq_vec, 'n', 1, width, nfft);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute significance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if num_shf ~= 0 
for i = 1:size(pacmat,1)
    for j = 1:size(pacmat,2)
        [h, p] = my_sig_test(pacmat(i,j), squeeze(shf_pacmat_final(:,i,j)), alpha);
        if h == 0
            pacmat(i,j) = 0;
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [ph_filt_signals, amp_filt_signals] = filt_signals(sig1,sig2,...
    Fs, ph_freq_vec, amp_freq_vec, measure, width, TFmethod)

total_num_dp = size(sig1,1);
num_trials = 1:size(sig1,2);

phfreq_low = min(ph_freq_vec);
phfreq_high = max(ph_freq_vec);
phfreq_bw = diff(ph_freq_vec(1:2));
ampfreq_low = min(amp_freq_vec);
ampfreq_high = max(amp_freq_vec);
ampfreq_bw = diff(amp_freq_vec(1:2));

xbins = ceil((phfreq_high - phfreq_low)/phfreq_bw);
ybins = ceil((ampfreq_high - ampfreq_low)/ampfreq_bw);


% Create structures to store filtered signals
ph_filt_signals = cell(1,xbins);
amp_filt_signals = cell(1,ybins);

for i2=1:xbins
    ph_filt_signals{1,i2} = zeros(total_num_dp,max(num_trials));
end


for i2=1:ybins
    amp_filt_signals{1,i2} = zeros(total_num_dp,max(num_trials));
end

% Filter and store filtered signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter sig2 for phase freq
for i2=1:xbins
    
    BandBounds(2) = phfreq_low+(i2*phfreq_bw);
    BandBounds(1) = BandBounds(2)-phfreq_bw;
     
    if strcmp(measure, 'esc')
        for i3=1:max(num_trials)
            ph_filt_signals{1,i2}(:,i3) = bp_vec(sig2(:,i3)', BandBounds, Fs, width, TFmethod)';
        end
    else
        for i3=1:max(num_trials)
            ph_filt_signals{1,i2}(:,i3) = ph_vec(sig2(:,i3)', BandBounds, Fs, width, TFmethod)';
        end
    end
    
       
end

% Filter sig1 for amplitude freq
for i2=1:ybins
    
    BandBounds(2) = ampfreq_low+(i2*ampfreq_bw);
    BandBounds(1) = BandBounds(2)-ampfreq_bw;
            
    for i3=1:max(num_trials)
        amp_filt_signals{1,i2}(:,i3) = amp_vec(sig1(:,i3)', BandBounds, Fs, width, TFmethod)';

    end
    
  
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
function shf_sig = shuffle_esc(sig)
% function shf_sig = shuffle_esc(sig, Fs)
%
% This function shuffles a signal using random insertion type shuffling. 
% The signal 'sig' is divided into a number of sections, equal to the
% number of seconds of data or 1000, which ever is smaller. These sections 
% are of random lengths, chosen from a uniform distribution. The sections 
% are then randomly re-ordered. 
% The input signal 'sig' must contain several seconds on data (>2) and must
% be a vector (i.e. number of trials = 1) or a cell array containing 
% multiple signals
%
% INPUTS:
% sig - input signal which is to be shuffled, passed as a column vector
% Fs - sampling frequency of 'sig', in Hz
%
% OUTPUTS:
% shf_sig - either a vector or a cell array depending on the input 'sig'.
% When used within find_pac_shf.m this output is of the same dimension as 
% filt_sig_mod and filt_sig_pac: number of cells - number of frequency bins
% and each cell element is a matrix(num_data_points, num_trials)
%
% Author: Rafal Bogacz, Angela Onslow, May 2010

sig_type = class(sig);

if iscell(sig)
    ybins = size(sig,1);
    xbins = size(sig,2);
    shf_sig_cell = cell(size(sig));   
else
    shf_sig = [];   
end
num_sec = 2;



switch sig_type
    
    case 'double'
        
        % Choose num_sec random 'cut' positions
        dpsplit = ceil(size(sig,1).*rand(num_sec,1));
        % Arrange these in ascending order
        dpsplit = sort (dpsplit);

        start(1)=1;
        start(2:num_sec)=dpsplit(1:num_sec-1);
        ending(1:num_sec-1)=dpsplit(1:num_sec-1)-1;
        ending(num_sec) =  size(sig,1);

        order = 1:num_sec;%randperm(num_sec);
    
        for c = 1:num_sec
        
            %shuffle the signal
            shf_sig = [shf_sig; sig(start(order(c)):ending(order(c)),:)];
        
        end


    case 'cell'
        for i = 1:ybins
            for j = 1:xbins
                
                current_sig = sig{i,j};
                
                % Choose num_sec random 'cut' positions
                dpsplit = ceil(size(sig{1,1},1).*rand(num_sec,1));
                % Arrange these in ascending order
                dpsplit = sort (dpsplit);

                start(1)=1;
                start(2:num_sec)=dpsplit(1:num_sec-1);
                ending(1:num_sec-1)=dpsplit(1:num_sec-1)-1;
                ending(num_sec) =  size(sig{1,1},1);

                order = randperm(num_sec);
    
                for c = 1:num_sec
        
                    %shuffle the signal
                    shf_sig_cell{i,j} = [shf_sig_cell{i,j}; current_sig(start(order(c)):ending(order(c)),:)];
        
                end
                
            end
        end
        
        shf_sig = shf_sig_cell;
end
end
function [escval, mxlag ]= esc_measure(ph_sig, amp_sig, avg)
% function escval = esc_measure(ph_sig, amp_sig, avg)
%
% Returns a value (or vector of values) for the ESC measure calculated
% between two signals. Signals may contain one of more trials. Multiple
% trials may be averaged so as to return one ESC value or a vector of
% ESC values calculated for each trial may be returned, depending on the 
% 'avg' argument. Signals should be passed as column vectors, multiple 
% trials stored as multiple columns.
%
% INPUTS:
% ph_sig - signal filtered for a lower, modulating frequency (e.g. theta
% band oscillations)
%
% amp_sig - signal filtered for a higher, modulated frequency (e.g. gamma
% band oscillations)
%
% avg - string, either 'y' or 'n', determines whether ESC values are
% averaged over trials or returned as a vector
%
% Author: Angela Onslow, May 2010
maxlags = 500;
escsum = 0;
if size(ph_sig, 2) ~= size(amp_sig, 2)
    sprintf('Error - Signals must have the same number of trials')
    return
end
num_trials = size(ph_sig, 2);

if strcmp(avg, 'y')
    
    %Average over trials using the Fisher transform
    for c = 1:num_trials
            
            %r = corrcoef(ph_sig(:,c), amp_sig(:,c));
            [r, lags] = xcorr(ph_sig(:,c),amp_sig(:,c),maxlags,'coeff');
            [val, ind] = max(r);%r(1,2);
            escsum = escsum + atanh(val);
            mxlag(c,1) = lags(ind);
            %escsum = escsum + atanh(r(1,2));
            
    end

    escsum = escsum/num_trials;
    escval = tanh(escsum);
    
else
    escval = zeros(num_trials,1);

    for i = 1:num_trials
        [r, lags] = xcorr(ph_sig(:,i), amp_sig(:,i),maxlags,'coeff');
        %escsum=max(xcorr(ph_sig(:,c),amp_sig(:,c),maxlags,'coeff'));
        [val, ind] = max(r);%r(1,2);
        escval(i,1) = val;
        mxlag(i,1) = lags(ind);
    end
end
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
ms = -maxlag:maxlag;
for t=1:length(ms)
    m=ms(t);
    z = corrcoef(X(m+maxlag+(1:K)),Y(maxlag+(1:K)-1));
    c(t) = z(1,2);
end
[val,ind]=max(c);
c=val;
lags = ms(ind);

end
