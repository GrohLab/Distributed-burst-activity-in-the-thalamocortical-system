function [digitalReward, digitalPunish, digitalLick] = intanADC(rhdFile)
% For better interpretation, convert analog noise signal in a digital
% trigger signal.
workingDir = fileparts(rhdFile);

if ~exist(fullfile(workingDir,'intanVariables.mat'),'file')
    fprintf('\nReading the intan variables...\n')
    read_Intan_RHD2000_file(rhdFile)
    fprintf('Loading the intan variables...\n')
    load(fullfile(workingDir,'intanVariables.mat'),'board_adc_data','board_dig_in_data','frequency_parameters')
else
    fprintf('\nLoading the intan variables...\n')
    load(fullfile(workingDir,'intanVariables.mat'),'board_adc_data','board_dig_in_data','frequency_parameters')
end

% Set sampling rate of ADC and digital channels to 10kHz
if frequency_parameters.board_adc_sample_rate ~= 10000
    convers = frequency_parameters.board_adc_sample_rate/10000;
    if floor(convers) ~= convers
        error("Sampling rate is not dividable by 10000")
    end
    board_adc_data = board_adc_data(:,1:convers:end);
    adc_sample_rate = 10000;
else
    adc_sample_rate = 10000;
end

if frequency_parameters.board_dig_in_sample_rate ~= 10000
    convers = frequency_parameters.board_dig_in_sample_rate/10000;
    if floor(convers) ~= convers
        error("Sampling rate is not dividable by 10000")
    end
    board_dig_in_data = board_dig_in_data(:,1:convers:end);
end

digitalPunish = zeros(2, size(board_adc_data,2));
digitalLick = zeros(2, size(board_adc_data,2));
digitalReward = board_dig_in_data;

timeFrame = 0.1*adc_sample_rate; % Defines 0.1 sec windows

for i = 1:2 % Channel 1 and 2 record the piezo lick signal for LP1 and LP2 respectively
    % Take the first derivative of lick signal
    fprintf('Converting piezo lick signals of channel %d...\n',i)
    dydx = gradient(board_adc_data(i,:)) ./ gradient(1:size(board_adc_data,2));
    dydx = lowpass(dydx,1000,10000); % Lowpass filter the signal to remove high amplitude noise from punishment
    threshold = (max(dydx)-median(dydx))/8;
    % Get the median value of analogous signal as the reference
    median_analog = median(dydx);
    % Turn analog signal into a digital trigger signal for better
    % interpretation
    idx = false(1,size(board_adc_data,2));
    
    % Checks whether a) the value itself is above median and b) previous
    % values crossed this value, since the lick signal is oscillating.
    for ii = 1:size(board_adc_data,2)
        if dydx(ii) > (median_analog+threshold) || dydx(ii) < (median_analog-threshold)
            idx(ii) = true;
        elseif ii > timeFrame && (any(dydx(ii-timeFrame:ii) > (median_analog+threshold)) ...
                || any(dydx(ii-timeFrame:ii) < (median_analog-threshold)))
            idx(ii) = true;
        end
    end
    % Sets values to 1, where lick events were registered.
    digitalLick(i,idx) = 1;
end



for i = 3:4 % Channel 3 and 4 record the noise amplifier signal for LP1 and LP2 respectively
    % Get the median value of analogous signal as the reference
    fprintf('Converting noise amplifier signals of channel %d...\n',i)
    median_analog = median(board_adc_data(i,:));
    % Turn analog signal into a digital trigger signal for better
    % interpretation
    idx = false(1,size(board_adc_data,2));
    
    % Preallocation for speed-up
    noise_data = board_adc_data(i,:);
    
    % Checks whether a) the value itself is above median and b) previous
    % values crossed this value, since the noise signal is oscillating.
    for ii = 1:size(board_adc_data,2)
        if noise_data(ii) > (median_analog+0.1) || noise_data(ii) < (median_analog-0.1)
            idx(ii) = true;
        elseif ii > timeFrame && (any(noise_data(ii-timeFrame:ii) > (median_analog+0.1)) ...
                || any(noise_data(ii-timeFrame:ii) < (median_analog-0.1)))
            idx(ii) = true;
        end
    end
    % Sets values to 1, where punishment occured.
    digitalPunish(i-2,idx) = 1;
end

save(fullfile(workingDir,'ADC_Data.mat'),'digitalReward', 'digitalPunish', 'digitalLick')

end