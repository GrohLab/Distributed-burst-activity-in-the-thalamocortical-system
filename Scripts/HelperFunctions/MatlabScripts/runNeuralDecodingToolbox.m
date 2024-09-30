%% Basic script for decoding behavioral stimuli based on neural data
function [RastFig, DecodeFig, num_cv_splits, num_units] = runNeuralDecodingToolbox(raster_data_dir,area,condition,spikeType, ...
    unitType,trialType,classifier_labels,classifierName,varargin)
options = inputParser;

checkArea = @(x) ischar(x) || iscell(x) || isempty(x);
checkCondition = @(x) any(validatestring(condition,{'Reward','Punishment','onlyFirstLick','Lick','WhiskerContact_left','WhiskerContact_right',...
    'WhiskerContact_onlyLeftFirst','WhiskerContact_onlyRightFirst'}));
checkSpikeType = @(x) any(validatestring(spikeType,{'allSpikes','burstSpikes','tonicSpikes'}));
checkUnitType = @(x) any(validatestring(unitType,{'allUnits','burstUnits','tonicUnits','bothUnits','noneUnits','touchModul'}));
checkTrialType = @(x) any(validatestring(trialType,{'allTrials','onlyGo','onlyNoGo','onlyNeutral','onlyNarrow','onlyWide','onlyIntermediate','onlyLick','onlyNoLick'}));
checkLabels = @(x) any(validatestring(classifier_labels,{'Go_NoGo_Neutral','Lick','Wide_Narrow_Intermediate'}));
checkClassifier = @(x) any(validatestring(classifierName,{'max_correlation_coefficient_CL','poisson_naive_bayes_CL','libsvm_CL'}));

addRequired(options,'raster_data_dir',@ischar);
addRequired(options,'area',checkArea);
addRequired(options,'condition',checkCondition);
addRequired(options,'spikeType',checkSpikeType);
addRequired(options,'unitType',checkUnitType);
addRequired(options,'trialType',checkTrialType);
addRequired(options,'classifier_labels',checkLabels);
addRequired(options,'classifierName',checkClassifier);
addParameter(options,'num_cv_splits',NaN,@isnumeric)
addParameter(options,'numUnits',NaN,@isnumeric)
addParameter(options,'numBootstrap',NaN,@isnumeric)
addParameter(options,'shuffle',true,@islogical)
addParameter(options,'trialCompare',{'narrow','wide'},@iscell)
addParameter(options,'outputFolder','',@ischar)

parse(options,raster_data_dir,area,condition,spikeType,unitType,...
    trialType,classifier_labels,classifierName, varargin{:})

raster_data_dir = options.Results.raster_data_dir;
area = options.Results.area;
condition = options.Results.condition;
spikeType = options.Results.spikeType;
unitType = options.Results.unitType;
trialType = options.Results.trialType;
classifier_labels = options.Results.classifier_labels;
classifierName = options.Results.classifierName;
num_cv_splits = options.Results.num_cv_splits;
numUnits = options.Results.numUnits;
numBootstrap = options.Results.numBootstrap;
shuffle = options.Results.shuffle;
trialCompare = options.Results.trialCompare;
outputFolder = options.Results.outputFolder;
if isempty(outputFolder)
    outputFolder = raster_data_dir;
end


assert(exist(raster_data_dir,'dir'),...
    'runNeuralDecoding:NodataDir',...
    'Your data path does not exist.')

if ischar(area)
    area = {area}; % Turn char values into single cell
end
for ar = 1:numel(area)
    assert(any(validatestring(area{ar},{'All','BC','VPM','POm','ZIv',...
        'BC-RS','VPM-RS','POm-RS','ZIv-RS','BC-FS','VPM-FS','POm-FS','ZIv-FS'})))
end

area_names = {'BC','VPM','POm','ZIv'};
area_colors = {'#377eb8','#4daf4a','#984ea3','#ff7f00'};

if isempty(area) || all(cellfun(@(x) any(strncmp(x, area,2)), area_names))
    area_name = 'allAreas';
elseif isscalar(area)
    area_name = area{1};
else
    area_name = strjoin(area,'&');
end

if isscalar(area)
    plotColor = area_colors{cellfun(@(x) contains(area, x), area_names)};
else
    plotColor = '#000000';
end
%%  Plot examplary raster and a PSTH for one neuron

FileInfo = cellfun(@(x) dir(fullfile(raster_data_dir, sprintf('%s_%s_%s*_%s_*_raster_data.mat',condition,spikeType,x,trialType))), area, 'UniformOutput', false);
FileInfo = vertcat(FileInfo{:});

if ~isempty(FileInfo)
    % load and plot random raster for sanity check
    rand_idx = randi(height(FileInfo),1,1);
    load(fullfile(FileInfo(rand_idx).folder,FileInfo(rand_idx).name),'raster_data','raster_site_info')
    trigger = raster_site_info.alignment_event_time-1;

    RastFig = figure;
    % plot the rasters
    subplot(1, 2, 1);
    imagesc(~raster_data);
    colormap gray;
    line([trigger trigger], get(gca, 'YLim'), 'color', [1 0 0]);  % put a line at the time when the stimulus was shown
    title(sprintf('Rasters of unit #%s (%s)',raster_site_info.unitId,raster_site_info.unitArea),'Interpreter','none');
    ylabel('Trials')
    xlabel('Time (ms)')


    % plot the PSTH
    subplot(1, 2, 2);
    hold on
    binSize = 10; % in msec
    nBins = size(raster_data,2)/binSize;
    raster_data_reshape = reshape(sum(raster_data,1),binSize,nBins);
    spikesPerTrial = sum(raster_data_reshape)/size(raster_data,1);
    bar((binSize:binSize:size(raster_data,2))-binSize/2,spikesPerTrial,'FaceColor',plotColor,'EdgeColor','none','BarWidth',1);
    line([trigger trigger], get(gca, 'YLim'), 'color', [1 0 0]);  % put a line at the time when the stimulus was shown
    annotation('textbox',[0.57 0.83 0.1 0.05],'String',{'bin size';[num2str(binSize) ' msec']},'FitBoxToText','on','EdgeColor','none','VerticalAlignment', 'bottom','FontSize',8)
    title('PSTH')
    ylabel('Spikes per trial per bin')
    xlabel('Time (ms)')
    axis([0 size(raster_data,2) get(gca, 'YLim')])

    set(gcf, 'position', [200   400   900   350])
else
    fprintf('\nNo analyzed raster data files in the respective folder and the picked area.\n')
    fprintf('Consider running the create_raster_data_files.m function first.\n')
    return
end

% initialize the random number generator
rng(sum(100*clock)); %#ok<CLOCK>
fprintf('\nInitializing the matlab random number generator to an aribrary clock value\n')

%%  Bin the data

if isnan(numUnits) && isnan(numBootstrap)
    save_prefix = fullfile(outputFolder,sprintf('Binned_data_%s_%s_%s_%s_%s',condition,spikeType,unitType,area_name,trialType));
else
    save_prefix = fullfile(outputFolder,sprintf('Binned_data_%s_%s_%s_%s_%s_bootstrap%03i',...
        condition,spikeType,unitType,area_name,trialType,numBootstrap));
end
bin_width = 100;
step_size = 25;
binned_data_file_name = [save_prefix '_' num2str(bin_width) 'ms_bins_' num2str(step_size) 'ms_sampled.mat'];

% Define which units to analyze
% For specific condition: fullfile(raster_data_dir,'condition*')
% For specific area: fullfile(raster_data_dir,'condition_area*')
DecodeFig = figure;
if ~exist(binned_data_file_name,'file')
    % This is where the responseType files are generated and stored
    stageFold = strsplit(raster_data_dir, '\');

    % Get the responseTypes.mat file (narrow&wide) to assess touch responses
    scriptFullPath = matlab.desktop.editor.getActiveFilename();
    try load(regexprep(scriptFullPath, 'Scripts.*', 'Scripts\\userDataPath.mat'), 'cohortPath');
    catch
        error('No userDataPath.mat file found. Run the userDataPath.m file first.')
    end
    responseDir = fullfile(cohortPath, 'Analysis-Figures\Burstiness-Scatter');
    responseFile = fullfile(responseDir,stageFold{end},sprintf('responseTypes_%s_%s.mat',condition,strjoin(trialCompare,'&')));

    switch unitType
        case 'allUnits'
            FileInfo = cellfun(@(x) dir(fullfile(raster_data_dir, sprintf('%s_%s_%s*_%s*',condition,spikeType,x,trialType))), area, 'UniformOutput', false);
            FileInfo = vertcat(FileInfo{:});
            % Input files as cell array
            binned_data_file_name = create_binned_data_from_raster_data(fullfile({FileInfo.folder},{FileInfo.name}), save_prefix, bin_width, step_size);
        otherwise
            % Load the responseTypes for assessing touch modulation
            if ~exist(responseFile,'file')
                % Flip the trial types, to see if that file does exist
                if exist(fullfile(responseDir,stageFold{end},sprintf('responseTypes_%s_%s.mat',condition,strjoin(flip(trialCompare),'&'))),'file')
                    responseFile = fullfile(responseDir,stageFold{end},sprintf('responseTypes_%s_%s.mat',condition,strjoin(flip(trialCompare),'&')));
                    load(responseFile,'responseTypes')
                else
                    fprintf('\nNo responseTypes.mat file found for trial types %s. Run the ApertureResponseTypes.m script first.\n', strjoin(trialCompare,' & '))
                    return
                end
            else
                load(responseFile,'responseTypes')
            end

            % Filter the responseTypes file for the given parameters
            switch unitType
                case 'noneUnits'
                    responseIdx = false(height(responseTypes),1);
                    for ar = 1:numel(area)
                        if contains(area{ar},'-FS')
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & ~responseTypes.PutExcitatory & ...
                                ~responseTypes.BurstSign & ~responseTypes.TonicSign);
                        elseif contains(area{ar},'-RS')
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & responseTypes.PutExcitatory & ...
                                ~responseTypes.BurstSign & ~responseTypes.TonicSign);
                        else
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & ...
                                ~responseTypes.BurstSign & ~responseTypes.TonicSign);
                        end
                    end
                    responseTypes = responseTypes(responseIdx,:);
                case 'burstUnits'
                    responseIdx = false(height(responseTypes),1);
                    for ar = 1:numel(area)
                        if contains(area{ar},'-FS')
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & ~responseTypes.PutExcitatory & ...
                                responseTypes.BurstSign);
                        elseif contains(area{ar},'-RS')
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & responseTypes.PutExcitatory & ...
                                responseTypes.BurstSign);
                        else
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & ...
                                responseTypes.BurstSign);
                        end
                    end
                    responseTypes = responseTypes(responseIdx,:);
                case 'tonicUnits'
                    responseIdx = false(height(responseTypes),1);
                    for ar = 1:numel(area)
                        if contains(area{ar},'-FS')
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & ~responseTypes.PutExcitatory & ...
                                responseTypes.TonicSign);
                        elseif contains(area{ar},'-RS')
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & responseTypes.PutExcitatory & ...
                                responseTypes.TonicSign);
                        else
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & ...
                                responseTypes.TonicSign);
                        end
                    end
                    responseTypes = responseTypes(responseIdx,:);
                case 'touchModul'
                    responseIdx = false(height(responseTypes),1);
                    for ar = 1:numel(area)
                        if contains(area{ar},'-FS')
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & ~responseTypes.PutExcitatory & ...
                                responseTypes.TouchModul & ~responseTypes.BurstSign);
                        elseif contains(area{ar},'-RS')
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & responseTypes.PutExcitatory & ...
                                responseTypes.TouchModul & ~responseTypes.BurstSign);
                        else
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & ...
                                responseTypes.TouchModul & ~responseTypes.BurstSign);
                        end
                    end
                    responseTypes = responseTypes(responseIdx,:);
                case 'bothUnits'
                    responseIdx = false(height(responseTypes),1);
                    for ar = 1:numel(area)
                        if contains(area{ar},'-FS')
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & ~responseTypes.PutExcitatory & ...
                                responseTypes.BurstSign & responseTypes.TonicSign);
                        elseif contains(area{ar},'-RS')
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & responseTypes.PutExcitatory & ...
                                responseTypes.BurstSign & responseTypes.TonicSign);
                        else
                            responseIdx = responseIdx | (cellfun(@(x) strncmp(x,area{ar},2), responseTypes.Area) & ...
                                responseTypes.BurstSign & responseTypes.TonicSign);
                        end
                    end
                    responseTypes = responseTypes(responseIdx,:);
            end

            % Select random number of units
            if ~isnan(numUnits)
                randInt = randperm(height(responseTypes), numUnits);
                responseTypes = responseTypes(sort(randInt),:);
            end

            % Access the respective raster data files
            rasterFiles = cell(height(responseTypes),1);
            for i = 1:height(responseTypes)
                tempFile = sprintf('%s_%s_%s*_%s_#%i_stage%i_ses%i_cl%s_raster_data.mat',condition,spikeType,responseTypes.Area{i},trialType,...
                    responseTypes.Animal{i},responseTypes.Stage{i},responseTypes.Session{i},responseTypes.UnitID{i});
                FileInfo = dir(fullfile(raster_data_dir,tempFile));
                rasterFiles{i} = FileInfo(1).name;
            end
            if ~isempty(rasterFiles)
                binned_data_file_name = create_binned_data_from_raster_data(fullfile(raster_data_dir,rasterFiles), save_prefix, bin_width, step_size);
            else
                textString = 'No units with given parameters.';
                % Display text, that there are no units to decode
                text(DecodeFig,0.5,0.5, textString, 'Units', 'normalized', 'FontSize', 8,...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                axis(f,'off')

                num_cv_splits = NaN;
                num_units = 0;
                return
            end

    end
end

%%  Calculate how many times each stimulus has been shown to each neuron

% To get reasonable results you usually need at least 5 repetitions of each condition
load(binned_data_file_name,'binned_labels','binned_data');  % load the binned data
num_units = numel(binned_data);

fprintf('\nTotal number of units: %d\n',num_units)
fprintf('A unit count below 5 can lead to inaccurate results.\n')

if isnan(num_cv_splits)
    num_sites_with_k_repeats = nan(1,60);
    for i = 0:59
        inds_of_sites_with_at_least_k_repeats = find_sites_with_k_label_repetitions(binned_labels.(classifier_labels), i);
        num_sites_with_k_repeats(i + 1) = length(inds_of_sites_with_at_least_k_repeats);
    end

    factor = 1;
    stopLoop = false;

    while ~stopLoop
        % This is the highest repetition, that represents a certain amount of units
        maxRepetitions = find(num_sites_with_k_repeats>=num_units*factor,1,'last');
        if maxRepetitions > 10
            % use no more than 10 cross-validation splits
            % (which means that 9 examples of each object are used for training and 1 example of each object is used for testing)
            num_cv_splits = 10;
            fprintf('\nChoosing a %d : 1 split ratio for training and test.\n',num_cv_splits-1)
            fprintf('%.f%% of all units fit this split count.\n\n',(num_sites_with_k_repeats(maxRepetitions)/num_units)*100)
            stopLoop = true;
        elseif maxRepetitions > 4
            num_cv_splits = maxRepetitions-1;
            fprintf('\nChoosing a %d : 1 split ratio for training and test.\n',num_cv_splits-1)
            fprintf('%.f%% of all units fit this split count.\n\n',(num_sites_with_k_repeats(maxRepetitions)/num_units)*100)
            stopLoop = true;
        elseif factor < 0.5
            fprintf('\nLess than 50%% of units have > 5 repetitions of the picked condition (%s)\n',classifier_labels)
            fprintf('As this might lead to inaccurate results, the script is aborted here.\n')
            return
        else
            factor = factor-0.1;
        end
    end
else
    inds_of_sites_with_at_least_k_repeats = find_sites_with_k_label_repetitions(binned_labels.(classifier_labels), num_cv_splits);
    num_sites_with_k_repeats = length(inds_of_sites_with_at_least_k_repeats);
    fprintf('\nSplit set by user to %d : 1 ratio for training and test.\n',num_cv_splits-1)
    fprintf('%.f%% of all units fit this split count.\n\n',(num_sites_with_k_repeats/num_units)*100)
end
%%  Create a datasource object and initialize the classifier

% %%% ATTENTION!! CHANGE
% num_cv_splits = 15;
% fprintf('\nATTENTION: Number of splits are overwritten here to %i splits\n',num_cv_splits)

% create the basic datasource object
ds = basic_DS(binned_data_file_name, classifier_labels,  num_cv_splits);
% Include only those units with enough trials
binned_fields = fieldnames(binned_labels);
idx = ismember(binned_fields,classifier_labels);
ds.sites_to_use = find_sites_with_k_label_repetitions(binned_labels.(binned_fields{idx}), num_cv_splits);

% other useful options:

% if using the Poison Naive Bayes classifier, load the data as spike counts by setting the load_data_as_spike_counts flag to 1
%ds = basic_DS(binned_data_file_name, specific_binned_labels_names,  num_cv_splits, 1);

% can have multiple repetitions of each label in each cross-validation split (which is a faster way to run the code that uses most of the data)
%ds.num_times_to_repeat_each_label_per_cv_split = 2;

% can do the decoding on a subset of labels
%ds.label_names_to_use =  {'kiwi', 'flower', 'guitar', 'hand'};


% Create a feature preprocessor object
% z-scoring normalizes each feature
the_feature_preprocessors{1} = zscore_normalize_FP;

% other useful options:

% can include a feature-selection features preprocessor to only use the top k most selective neurons
% fp = select_or_exclude_top_k_features_FP;
% fp.num_features_to_use = 25;   % use only the 25 most selective neurons as determined by a univariate one-way ANOVA
% the_feature_preprocessors{2} = fp;


% Create a classifier object
% note that for the poisson naive bayes classifier the data needs to be loaded as spike counts
the_classifier = eval(classifierName);

if isequal(classifierName,'libsvm_CL')
    %     the_classifier.kernel = 'polynomial';
    %     the_classifier.poly_degree = 3;
    the_classifier.kernel = 'gaussian';
    the_classifier.gaussian_gamma = 3;
end

% Create the cross-validator
the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
the_cross_validator.num_resample_runs = 10;  % usually more than 2 resample runs are used to get more accurate results, but to save time we are using a small number here
the_cross_validator.display_progress.zero_one_loss = 0;

% other useful options:
% can greatly speed up the run-time of the analysis by not creating a full TCT matrix (i.e., only trainging and testing the classifier on the same time bin)
% the_cross_validator.test_only_at_training_times = 1;

%%  Run decoding analysis and save results

% if calling the code from a script, one can log the code so that one can recreate the results
%log_code_obj = log_code_object;
%log_code_obj.log_current_file;
% if logged the code that was run using a log_code_object, save the code
%LOGGED_CODE = log_code_obj.return_logged_code_structure;
%save(save_file_name, '-v7.3', 'DECODING_RESULTS', 'LOGGED_CODE');

if isnan(numUnits) && isnan(numBootstrap)
    save_file_name = fullfile(outputFolder,[sprintf('Binned_data_results_%s_%s_%s_%s_%s_%s_%isplits',...
        condition,spikeType,unitType,area_name,trialType,classifierName,num_cv_splits),'.mat']);
else
    save_file_name = fullfile(outputFolder,[sprintf('Binned_data_results_%s_%s_%s_%s_%s_%s_%isplits_bootstrap%03i',...
        condition,spikeType,unitType,area_name,trialType,classifierName,num_cv_splits,numBootstrap),'.mat']);
end

if ~exist(save_file_name,'file')
    % run the decoding analysis
    DECODING_RESULTS = the_cross_validator.run_cv_decoding;
    save(save_file_name, 'DECODING_RESULTS');
end

%% Run shuffled data for calculating p-values
if shuffle
    shuffle_dir_name = fullfile(outputFolder,sprintf('Shuffle_results_%s_%s_%s_%s_%s_%s_%isplits',...
        condition,classifierName,spikeType,unitType,area_name,trialType,num_cv_splits));
    if exist(shuffle_dir_name,'dir') == 0
        mkdir(shuffle_dir_name)
        fprintf('\nRunning shuffled decoding for statistical comparison.\n')
        for shuff_num = 1:5
            fprintf('Shuffle iteration: %d/5\n',shuff_num)
            ds_shuff = basic_DS(binned_data_file_name, classifier_labels,  num_cv_splits);
            ds_shuff.sites_to_use = find_sites_with_k_label_repetitions(binned_labels.(binned_fields{idx}), num_cv_splits);
            ds_shuff.randomly_shuffle_labels_before_running=1;

            the_cross_validator = standard_resample_CV(ds_shuff, the_classifier, the_feature_preprocessors);
            the_cross_validator.num_resample_runs = 10;

            % suppress displays
            the_cross_validator.display_progress.zero_one_loss = 0;
            the_cross_validator.display_progress.resample_run_time = 0;

            DECODING_RESULTS = the_cross_validator.run_cv_decoding;

            save(fullfile(shuffle_dir_name,sprintf('ShuffRun_%02d.mat',shuff_num)),...
                'DECODING_RESULTS');
        end
    end

%%  Plot the basic results
load(save_file_name, 'DECODING_RESULTS')

% which results should be plotted (only have one result to plot here)
result_names{1} = save_file_name;
pval_dir_name{1} = [shuffle_dir_name,'\'];

% create an object to plot the results
plot_obj = plot_standard_results_object(result_names);

plot_obj.the_colors = {hex2rgb(plotColor)}; % must be cell and rgb format
plot_obj.errorbar_file_names = result_names;
plot_obj.errorbar_transparency_level = .15;
plot_obj.errorbar_edge_transparency_level = 0.05;

numSplits = DECODING_RESULTS.DS_PARAMETERS.num_cv_splits;
plot_obj.errorbar_stdev_multiplication_factor = 1/sqrt(numSplits); %to make it SEM

plot_obj.p_values = pval_dir_name;
plot_obj.collapse_all_times_when_estimating_pvals = 1;
plot_obj.p_value_alpha_level = 1e-10;
plot_obj.add_pvalue_latency_to_legends_alignment = 0;

% plot_obj.chance_level = 0.5;

% put a line at the time when the stimulus was shown
plot_obj.significant_event_times = 0;
plot_obj.the_axis = [-inf inf -inf 100];

% optional argument, can plot different types of results
% plot_obj.result_type_to_plot = 2;  % for example, setting this to 2 plots the normalized rank results

% plot the results as usual
plot_obj.plot_results;

ylabel('Accuracy [%]')
xlabel('Time [msec]')


switch classifier_labels
    case 'Go_NoGo_Neutral'
        title({sprintf('Trial type prediction accuracy based on %s',condition);sprintf('%s - %s',area_name,spikeType)},'Interpreter','none')
    case 'Lick'
        title({sprintf('Lick prediction accuracy based on %s',condition);sprintf('%s - %s',area_name,spikeType)},'Interpreter','none')
    case 'Wide_Narrow_Intermediate'
        title({sprintf('Aperture prediction accuracy based on %s',condition);sprintf('%s - %s',area_name,spikeType)},'Interpreter','none')
end

savefig(fullfile(raster_data_dir,sprintf('plot_%s_%s_%s_%s_%s_%s_%s_%dsplits',condition,classifier_labels,spikeType,unitType,area_name,trialType,classifierName,numSplits)))

end

end