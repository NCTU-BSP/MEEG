function varargout = process_SAE( varargin )
% PROCESS_TTEST: Student''s t-test: Compare means between conditions (across trials or across sujects).

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2015 University of Southern California & McGill University
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
% Authors: Francois Tadel, Dimitrios Pantazis, 2008-2013
macro_methodcall;
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment = 'Stacked Autoendcoder';
    sProcess.FileTag = '';
    sProcess.Category    = 'Stat2';
    sProcess.SubGroup    = 'Test';
    sProcess.Index       = 1000;
%     % Definition of the input accepted by this process
     sProcess.InputTypes  = {'data', 'results', 'timefreq', 'matrix'};
     sProcess.OutputTypes = {'data', 'results', 'timefreq', 'matrix'};
     sProcess.nInputs     = 2;
     sProcess.nMinFiles   = 2;
     % Hidden Layer1 
    sProcess.options.label1.Comment = '<HTML><BR><U><B>Hidden Layer1</B></U>:';
    sProcess.options.label1.Type    = 'label';
    sProcess.options.hid1_neuron.Comment = 'Number of hidden layer 1 neurons: ';
    sProcess.options.hid1_neuron.Type = 'value';
    sProcess.options.hid1_neuron.Value = {100,'',0};
    sProcess.options.hid1_learningRate.Comment = 'hidden layer 1 learning rate: ';
    sProcess.options.hid1_learningRate.Type = 'value';
    sProcess.options.hid1_learningRate.Value = {0.01,'',4};
    
    % Hidden Layer2
    sProcess.options.label2.Comment = '<HTML><BR><U><B>Hidden Layer2</B></U>:';
    sProcess.options.label2.Type    = 'label';
    sProcess.options.hid2_neuron.Comment = 'Number of hidden layer 2 neurons: ';
    sProcess.options.hid2_neuron.Type = 'value';
    sProcess.options.hid2_neuron.Value = {100,'',0};   
    sProcess.options.hid2_learningRate.Comment = 'hidden layer 2 learning rate: ';
    sProcess.options.hid2_learningRate.Type = 'value';
    sProcess.options.hid2_learningRate.Value = {0.01,'',4};
    
    
    sProcess.options.label3.Comment = '<HTML><BR><U><B>Model Parameter</B></U>:';
    sProcess.options.label3.Type    = 'label';
    % number of epoch
    sProcess.options.niter.Comment = 'Number of iteraitions';
    sProcess.options.niter.Type = 'value';
    sProcess.options.niter.Value = {100,'',0};
    
    % number of batchsize
    sProcess.options.batchsize.Comment = 'Number of batch size';
    sProcess.options.batchsize.Type = 'value';
    sProcess.options.batchsize.Value = {[],'',0};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
 Comment = 'SAE MODEL';
end


%% ===== RUN =====
function sOutput = Run(sProcess, sInputsA, sInputsB) %#ok<DEFNU>

sInputsA(1)

% process inputA
fiA_ln = length(sInputsA);
data_A = [];

 for i = 1:fiA_ln 
    datamat = in_bst(sInputsA(i).FileName);
    if strcmp(sInputsA(i).FileType,'timefreq') == 1 
        size_data = size(squeeze(datamat.TF));
        single_data = reshape(squeeze(datamat.TF),1,size_data(1)*size_data(2));
        data_A = vertcat(data_A,single_data);      
    else
        size_data = size(datamat.F);
        single_data = reshape(datamat.F,1,size_data(1)*size_data(2));
        data_A = vertcat(data_A,single_data);
    end
 end
 size_data
% process inputB
fiB_ln = length(sInputsB);
data_B = [];

 for i = 1:fiB_ln 
    datamat = in_bst(sInputsB(i).FileName);
    if strcmp(sInputsA(i).FileType,'timefreq') == 1 
        size_data = size(squeeze(datamat.TF));
        single_data = reshape(squeeze(datamat.TF),1,size_data(1)*size_data(2));
        data_B = vertcat(data_B,single_data);      
    else
        size_data = size(datamat.F);
        single_data = reshape(datamat.F,1,size_data(1)*size_data(2));
        data_B = vertcat(data_B,single_data);
    end
 end
 
 allA_ln = size(data_A);
 allB_ln = size(data_B);
 
 % make label
 label = zeros(allA_ln(1)+allB_ln(1),2);
 label(1:allA_ln(1),1)=1;
 label(allA_ln(1)+1:end,2)=1;

 alldata = [data_A;data_B];
 all_ln = size(alldata);
 %%----------------SAE model-----------------------
 rand('state',0)
 sae = saesetup([all_ln(2) sProcess.options.hid1_neuron.Value{1} sProcess.options.hid2_neuron.Value{1}]);
 sae.ae{1}.activation_function       = 'sigm';
 sae.ae{2}.activation_function       = 'sigm';
 sae.ae{1}.learningRate              =  sProcess.options.hid1_learningRate.Value{1};
 sae.ae{2}.learningRate              =  sProcess.options.hid2_learningRate.Value{1};
 
opts.numepochs = sProcess.options.niter.Value{1};
opts.batchsize = sProcess.options.batchsize.Value{1};
sae = saetrain(sae, alldata, opts);

nn = nnsetup([all_ln(2) sProcess.options.hid1_neuron.Value{1} sProcess.options.hid2_neuron.Value{1} 2]);
nn.activation_function              = 'sigm';
nn.learningRate                     = sProcess.options.hid1_learningRate.Value{1};
nn.W{1} = sae.ae{1}.W{1};
nn.W{2} = sae.ae{2}.W{1};
nn = nntrain(nn, alldata, label, opts);

save('sae_model','nn');
sOutput = [];
end

%%%%%%%%%%%  SAE  functions form deep learn toolbox, 2012, Rasmus %%%%%%%%%%%%%%%

function sae = saesetup(size)
    for u = 2 : numel(size)
        sae.ae{u-1} = nnsetup([size(u-1) size(u) size(u-1)]);
    end
end

function sae = saetrain(sae, x, opts)
    for i = 1 : numel(sae.ae);
        disp(['Training AE ' num2str(i) '/' num2str(numel(sae.ae))]);
        sae.ae{i} = nntrain(sae.ae{i}, x, x, opts);
        t = nnff(sae.ae{i}, x, x);
        x = t.a{2};
        %remove bias term
        x = x(:,2:end);
    end
end

function nn = nnsetup(architecture)
%NNSETUP creates a Feedforward Backpropagate Neural Network
% nn = nnsetup(architecture) returns an neural network structure with n=numel(architecture)
% layers, architecture being a n x 1 vector of layer sizes e.g. [784 100 10]

    nn.size   = architecture;
    nn.n      = numel(nn.size);
    
    nn.activation_function              = 'sigm';   %  Activation functions of hidden layers: 'sigm' (sigmoid) or 'tanh_opt' (optimal tanh).
    nn.learningRate                     = 2;            %  learning rate Note: typically needs to be lower when using 'sigm' activation function and non-normalized inputs.
    nn.momentum                         = 0.5;          %  Momentum
    nn.scaling_learningRate             = 1;            %  Scaling factor for the learning rate (each epoch)
    nn.weightPenaltyL2                  = 0;            %  L2 regularization
    nn.nonSparsityPenalty               = 0;            %  Non sparsity penalty
    nn.sparsityTarget                   = 0.05;         %  Sparsity target
    nn.inputZeroMaskedFraction          = 0;            %  Used for Denoising AutoEncoders
    nn.dropoutFraction                  = 0;            %  Dropout level (http://www.cs.toronto.edu/~hinton/absps/dropout.pdf)
    nn.testing                          = 0;            %  Internal variable. nntest sets this to one.
    nn.output                           = 'sigm';       %  output unit 'sigm' (=logistic), 'softmax' and 'linear'

    for i = 2 : nn.n   
        % weights and weight momentum
        nn.W{i - 1} = (rand(nn.size(i), nn.size(i - 1)+1) - 0.5) * 2 * 4 * sqrt(6 / (nn.size(i) + nn.size(i - 1)));
        nn.vW{i - 1} = zeros(size(nn.W{i - 1}));
        
        % average activations (for use with sparsity)
        nn.p{i}     = zeros(1, nn.size(i));  

    end
end

function [nn, L]  = nntrain(nn, train_x, train_y, opts, val_x, val_y)
%NNTRAIN trains a neural net
% [nn, L] = nnff(nn, x, y, opts) trains the neural network nn with input x and
% output y for opts.numepochs epochs, with minibatches of size
% opts.batchsize. Returns a neural network nn with updated activations,
% errors, weights and biases, (nn.a, nn.e, nn.W, nn.b) and L, the sum
% squared error for each training minibatch.

assert(isfloat(train_x), 'train_x must be a float');
assert(nargin == 4 || nargin == 6,'number ofinput arguments must be 4 or 6')

loss.train.e               = [];
loss.train.e_frac          = [];
loss.val.e                 = [];
loss.val.e_frac            = [];
opts.validation = 0;
if nargin == 6
    opts.validation = 1;
end

fhandle = [];
if isfield(opts,'plot') && opts.plot == 1
    fhandle = figure();
end

m = size(train_x, 1);

batchsize = opts.batchsize;
numepochs = opts.numepochs;

numbatches = m / batchsize;

assert(rem(numbatches, 1) == 0, 'numbatches must be a integer');

L = zeros(numepochs*numbatches,1);
%L = gdouble(L);
n = 1;

for i = 1 : numepochs
    tic;
    
    kk = randperm(m);
    for l = 1 : numbatches
        batch_x = train_x(kk((l - 1) * batchsize + 1 : l * batchsize), :);
        
        %Add noise to input (for use in denoising autoencoder)
        if(nn.inputZeroMaskedFraction ~= 0)
            batch_x = batch_x.*(rand(size(batch_x))>nn.inputZeroMaskedFraction);
%             batch_x = batch_x.*(grand(size(batch_x))>nn.inputZeroMaskedFraction);  %GPU
        end
        
        batch_y = train_y(kk((l - 1) * batchsize + 1 : l * batchsize), :);
        
        nn = nnff(nn, batch_x, batch_y);
        nn = nnbp(nn);
        nn = nnapplygrads(nn);
        
%         L(n) = gdouble(nn.L); %GPU
        L(n) = nn.L;
        n = n + 1;
    end
    
    t = toc;
    
     if opts.validation == 1
        loss = nneval(nn, loss, train_x, train_y, val_x, val_y);
        str_perf = sprintf('; Full-batch train mse = %f, val mse = %f', loss.train.e(end), loss.val.e(end));
     else
        loss = nneval(nn, loss, train_x, train_y);
        str_perf = sprintf('; Full-batch train err = %f', loss.train.e(end));
     end

      disp(['epoch ' num2str(i) '/' num2str(opts.numepochs) '. Took ' num2str(t) ' seconds' '. Mini-batch mean squared error on training set is ' num2str(mean(L((n-numbatches):(n-1)))) str_perf]);
      nn.learningRate = nn.learningRate * nn.scaling_learningRate;

end
end

function nn = nnapplygrads(nn)
%NNAPPLYGRADS updates weights and biases with calculated gradients
% nn = nnapplygrads(nn) returns an neural network structure with updated
% weights and biases
    
    for i = 1 : (nn.n - 1)
        if(nn.weightPenaltyL2>0)
            %dW = nn.dW{i} + gdouble(nn.weightPenaltyL2) * [gzeros(size(nn.W{i},1),1) nn.W{i}(:,2:end)]; %GPU
            dW = nn.dW{i} + nn.weightPenaltyL2 * [zeros(size(nn.W{i},1),1) nn.W{i}(:,2:end)];
        else
            dW = nn.dW{i};
        end
        
        %dW = gdouble(nn.learningRate) * dW;
        dW = nn.learningRate * dW;
        if(nn.momentum>0)
            %nn.vW{i} = gdouble(nn.momentum)*nn.vW{i} + dW;
            nn.vW{i} = nn.momentum*nn.vW{i} + dW;
            dW = nn.vW{i};
        end
            
        nn.W{i} = nn.W{i} - dW;
    end
end

function nn = nnbp(nn)
%NNBP performs backpropagation
% nn = nnbp(nn) returns an neural network structure with updated weights 
    
    n = nn.n;
    sparsityError = 0;
    switch nn.output
        case 'sigm'
            d{n} = - nn.e .* (nn.a{n} .* (1 - nn.a{n}));
        case {'softmax','linear'}
            d{n} = - nn.e;
    end
    for i = (n - 1) : -1 : 2
        % Derivative of the activation function
        switch nn.activation_function 
            case 'sigm'
                d_act = nn.a{i} .* (1 - nn.a{i});
            case 'tanh_opt'
                d_act = 1.7159 * 2/3 * (1 - 1/(1.7159)^2 * nn.a{i}.^2);
        end
        
        if(nn.nonSparsityPenalty>0)
            pi = repmat(nn.p{i}, size(nn.a{i}, 1), 1);
            sparsityError = [zeros(size(nn.a{i},1),1) nn.nonSparsityPenalty * (-nn.sparsityTarget ./ pi + (1 - nn.sparsityTarget) ./ (1 - pi))];
        end
        
        % Backpropagate first derivatives
        if i+1==n % in this case in d{n} there is not the bias term to be removed             
            d{i} = (d{i + 1} * nn.W{i} + sparsityError) .* d_act; % Bishop (5.56)
        else % in this case in d{i} the bias term has to be removed
            d{i} = (d{i + 1}(:,2:end) * nn.W{i} + sparsityError) .* d_act;
        end
        
        if(nn.dropoutFraction>0)
            d{i} = d{i} .* [ones(size(d{i},1),1) nn.dropOutMask{i}];
        end

    end

    for i = 1 : (n - 1)
        if i+1==n
            nn.dW{i} = (d{i + 1}' * nn.a{i}) / size(d{i + 1}, 1);
        else
            nn.dW{i} = (d{i + 1}(:,2:end)' * nn.a{i}) / size(d{i + 1}, 1);      
        end
    end
 end

 function nn = nnff(nn, x, y)
%NNFF performs a feedforward pass
% nn = nnff(nn, x, y) returns an neural network structure with updated
% layer activations, error and loss (nn.a, nn.e and nn.L)

    n = nn.n;
    m = size(x, 1);
    
    x = [ones(m,1) x];
%     nn.a{1} = gdouble(x); %GPU
    nn.a{1} = x;
    %feedforward pass
    for i = 2 : n-1
        switch nn.activation_function 
            case 'sigm'
                % Calculate the unit's outputs (including the bias term)
                nn.a{i} = sigm(nn.a{i - 1} * nn.W{i - 1}');
            case 'tanh_opt'
                nn.a{i} = tanh_opt(nn.a{i - 1} * nn.W{i - 1}');
        end
        
        %dropout
        if(nn.dropoutFraction > 0)
            if(nn.testing)
                nn.a{i} = nn.a{i}.*(1 - nn.dropoutFraction);
            else
                nn.dropOutMask{i} = (rand(size(nn.a{i}))>nn.dropoutFraction);
                nn.a{i} = nn.a{i}.*nn.dropOutMask{i};
            end
        end
        
        %calculate running exponential activations for use with sparsity
        if(nn.nonSparsityPenalty>0)
            nn.p{i} = 0.99 * nn.p{i} + 0.01 * mean(nn.a{i}, 1);
        end
        
        %Add the bias term
        nn.a{i} = [ones(m,1) nn.a{i}];
    end
    switch nn.output 
        case 'sigm'
            nn.a{n} = sigm(nn.a{n - 1} * nn.W{n - 1}');
        case 'linear'
            nn.a{n} = nn.a{n - 1} * nn.W{n - 1}';
        case 'softmax'
            nn.a{n} = nn.a{n - 1} * nn.W{n - 1}';
            nn.a{n} = exp(bsxfun(@minus, nn.a{n}, max(nn.a{n},[],2)));
            nn.a{n} = bsxfun(@rdivide, nn.a{n}, sum(nn.a{n}, 2)); 
    end

    %error and loss
    nn.e = (y - nn.a{n});
    
    switch nn.output
        case {'sigm', 'linear'}
            nn.L = 1/2 * sum(sum(nn.e .^ 2)) / m; 
        case 'softmax'
            nn.L = -sum(sum(y .* log(nn.a{n}))) / m;
    end
 end

 function [loss] = nneval(nn, loss, train_x, train_y, val_x, val_y)
%NNEVAL evaluates performance of neural network
% Returns a updated loss struct
assert(nargin == 4 || nargin == 6, 'Wrong number of arguments');

nn.testing = 1;
% training performance
nn                    = nnff(nn, train_x, train_y);
loss.train.e(end + 1) = nn.L;

% validation performance
if nargin == 6
    nn                    = nnff(nn, val_x, val_y);
    loss.val.e(end + 1)   = nn.L;
end
nn.testing = 0;
%calc misclassification rate if softmax
if strcmp(nn.output,'softmax')
    [er_train, dummy]               = nntest(nn, train_x, train_y);
    loss.train.e_frac(end+1)    = er_train;
    
    if nargin == 6
        [er_val, dummy]             = nntest(nn, val_x, val_y);
        loss.val.e_frac(end+1)  = er_val;
    end
end

 end

function [er, bad] = nntest(nn, x, y)
    labels = nnpredict(nn, x);
    [dummy, expected] = max(y,[],2);
    bad = find(labels ~= expected);    
    er = numel(bad) / size(x, 1);
end

function labels = nnpredict(nn, x)
    nn.testing = 1;
    nn = nnff(nn, x, zeros(size(x,1), nn.size(end)));
    nn.testing = 0;
    
    [dummy, i] = max(nn.a{end},[],2);
    labels = i;
end

function X = sigm(P)
    X = 1./(1+exp(-P));
end