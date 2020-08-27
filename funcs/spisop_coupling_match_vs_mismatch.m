%% Script to get SpiSOP coupling output to topoplot format

% Create one structure per participant containing the average
% amount of co-occurences, the mean delay and SD delay per channel.
% In addition, filter out the multiple spindle occurences per target window
% by removing any duplicates

clear; clear all
addpath(genpath('...\spisop\output\...\'))

channel_info = {'EEG_-C3/A2'; 'EEG_-C4/A1'};      
     

%% Matching events
%%% Read in the SpiSOP coupling output data
subject_numbers = readtable('events_cooccurrance_matching_all_recent.csv');
%%% find unique subject numbers
subject_numbers = unique(subject_numbers.event_files_number);

for sub_nr = 1:size(subject_numbers,1)
    tic
    %%% Read in matching events for fast spindles
    input_name = sprintf('events_cooccurrance_matching_events_file_num_%d.csv', subject_numbers(sub_nr));
    soe_coupling = readtable(input_name);
 
    %%% Tag the duplicates with a 0
    for channel_cntr = 1:height(soe_coupling)-1
        if strcmp(soe_coupling.test_channel(channel_cntr),(soe_coupling.test_channel(channel_cntr+1))) == 1 %% match channel name
            if soe_coupling.test_seconds_begin(channel_cntr) == soe_coupling.test_seconds_begin(channel_cntr+1) %% match spindle onset
            soe_coupling.test_seconds_begin(channel_cntr+1) = 0; % Tag duplicate
            soe_coupling.test_channel(channel_cntr+1) = {0}; % Tag duplicate (redundant)
            end
        end
    end

    %%% Remove the tagged duplicates
    duplicates_index = soe_coupling.test_seconds_begin == 0;
    soe_coupling(duplicates_index,:) = [];

 %%% Extract amount of co-occurences per subject
    %%% Create a subject specific channel_info vector
    channel_info_subj = channel_info;

    %%% Extract channel names and position
    [C,ia,ic] = unique(soe_coupling.test_channel, 'stable');
    %%% Extract counts per channel
    a_counts = accumarray(ic,1);
    %%% Extract position of channel in vector to sort in correct order
    [tf,idx] = ismember(channel_info_subj(1,:),C);
    sum_coupling = num2cell(a_counts(idx)');
    %%% In case a cell is empty, change to 0
    empty_cells_index = cellfun('isempty', sum_coupling);
    sum_coupling(empty_cells_index) = {0};
    channel_info_subj(2,:) = sum_coupling;
    
 %%% Extract mean delay
    for k = 1:size(channel_info_subj,2)
        mask = strcmp(soe_coupling.test_channel, channel_info_subj(1,k));
        mean_delay_store = mean(soe_coupling.test_seconds_begin(mask) - soe_coupling.target_seconds_begin(mask));
        channel_info_subj{3,k}= mean_delay_store;
    end
    mean_delay =  channel_info_subj(3,:);
    
 %%% Extract sd delay
    for k = 1:size(channel_info_subj,2)
        mask = strcmp(soe_coupling.test_channel, channel_info_subj(1,k));
        sd_delay_store = std(soe_coupling.test_seconds_begin(mask) - soe_coupling.target_seconds_begin(mask));
        channel_info_subj{4,k}= sd_delay_store;
    end
    sd_delay =  channel_info_subj(4,:);
    
 %%% Extract all the spindle parameters
   spin_duration  = get_mean(soe_coupling.test_duration_seconds, soe_coupling.test_channel, channel_info_subj(1,:));
   spin_frequency = get_mean(soe_coupling.test_frequency_by_mean_pk_trgh_cnt_per_dur, soe_coupling.test_channel, channel_info_subj(1,:));
   spin_amplitude = get_mean(soe_coupling.test_amplitude_peak2trough_max, soe_coupling.test_channel, channel_info_subj(1,:));
    
 %%% Extract all the SO parameters
   so_amplitude = get_mean(soe_coupling.target_amplitude_peak2trough_max, soe_coupling.test_channel, channel_info_subj(1,:));
   so_slope     = get_mean(soe_coupling.target_slope_trough_to_zeroxing_potential_per_second, soe_coupling.test_channel, channel_info_subj(1,:));
   so_duration  = get_mean(soe_coupling.target_duration_seconds, soe_coupling.test_channel, channel_info_subj(1,:));

%% Mismatched events

%%% Read in mismatching events
    input_name = sprintf('events_cooccurrance_mismatching_events_file_num_%d.csv', subject_numbers(sub_nr));
    soe_mismatch = readtable(input_name);
    
     
%%% Create a subject specific channel_info vector
    channel_info_subj_mismatch = channel_info;

  %%% Extract amount of mismatches per subject
    %%% Extract channel names and position
    [C,ia,ic] = unique(soe_mismatch.test_channel, 'stable');
    %%% Extract counts per channel
    a_counts = accumarray(ic,1);
    %%% Extract position of channel in vector to sort in correct order
    [tf,idx] = ismember(channel_info_subj_mismatch(1,:),C);
    sum_mismatches = num2cell(a_counts(idx)');
    %%% In case a cell is empty, change to 0
    empty_cells_index = cellfun('isempty', sum_mismatches);
    sum_mismatches(empty_cells_index) = {0};
       
   % Extract all the spindle parameters
   spin_duration_mismatch  = get_mean(soe_mismatch.test_duration_seconds, soe_mismatch.test_channel, channel_info_subj(1,:));
   spin_frequency_mismatch = get_mean(soe_mismatch.test_frequency_by_mean_pk_trgh_cnt_per_dur, soe_mismatch.test_channel, channel_info_subj(1,:));
   spin_amplitude_mismatch = get_mean(soe_mismatch.test_amplitude_peak2trough_max, soe_mismatch.test_channel, channel_info_subj(1,:));

 % -------------------------------------------------------------------------------------------------------------
     %%% Collect mismatching SO data  
     %%% Read in so events
      so_input_name = sprintf('nonsens_run_1_nonsens_so_events_datanum_%d.csv', subject_numbers(sub_nr));
      so_all = readtable(so_input_name);
  
    %%% Tag the duplicates with a 0
    so_mismatch_mask = true(size(so_all,1),1);
      for so_coup_counter = 1:height(soe_coupling)
          idx = find(so_all.seconds_begin==soe_coupling.target_seconds_begin(so_coup_counter) &...
              strcmp(so_all.channel, soe_coupling.target_channel{so_coup_counter}));         
          so_mismatch_mask(idx) = 0;
      end
    so_mismatch = so_all(so_mismatch_mask, :);
    %%% Remove the tagged duplicates
    match_index = so_mismatch.datasetnum == 0;
    so_mismatch(match_index,:) = [];

  %%% Extract amount of SO mismatches per subject
    %%% Extract channel names and position
    [C,ia,ic] = unique(so_mismatch.channel, 'stable');
    %%% Extract counts per channel
    a_counts = accumarray(ic,1);
    %%% Extract position of channel in vector to sort in correct order
    [tf,idx] = ismember(channel_info_subj_mismatch(1,:),C);
    so_mismatch_amount = num2cell(a_counts(idx)');
    %%% In case a cell is empty, change to 0
    empty_cells_index = cellfun('isempty', so_mismatch_amount);
    so_mismatch_amount(empty_cells_index) = {0};
    
  %%% Extract amount of SO (total) per subject
    %%% Extract channel names and position
    [C,ia,ic] = unique(so_all.channel, 'stable');
    %%% Extract counts per channel
    a_counts = accumarray(ic,1);
    %%% Extract position of channel in vector to sort in correct order
    [tf,idx] = ismember(channel_info_subj(1,:),C);
    so_all_amount = num2cell(a_counts(idx)');
    %%% In case a cell is empty, change to 0
    empty_cells_index = cellfun('isempty', so_all_amount);
    so_all_amount(empty_cells_index) = {0};
    
 %%% Extract all the SO parameters
   so_amplitude_mismatch = get_mean(so_mismatch.amplitude_peak2trough_max, so_mismatch.channel, channel_info_subj(1,:));
   so_slope_mismatch     = get_mean(so_mismatch.slope_trough_to_zeroxing_potential_per_second, so_mismatch.channel, channel_info_subj(1,:));
   so_duration_mismatch  = get_mean(so_mismatch.duration_seconds, so_mismatch.channel, channel_info_subj(1,:));
    
  %%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
    %%% Create structures
    coupling.match.sum                   = sum_coupling;
    coupling.mismatch.sum                = sum_mismatches;
    coupling.match.mean_delay            = mean_delay;
    coupling.match.sd_delay              = sd_delay;
    
    coupling.match.spindle.duration      = spin_duration;
    coupling.match.spindle.frequency     = spin_frequency;
    coupling.match.spindle.amplitude     = spin_amplitude;
    coupling.mismatch.spindle.duration   = spin_duration_mismatch;
    coupling.mismatch.spindle.frequency  = spin_frequency_mismatch;
    coupling.mismatch.spindle.amplitude  = spin_amplitude_mismatch;    

    coupling.match.so.amplitude          = so_amplitude;
    coupling.match.so.slope              = so_slope;
    coupling.match.so.duration           = so_duration;
    coupling.mismatch.so.amplitude       = so_amplitude_mismatch;
    coupling.mismatch.so.slope           = so_slope_mismatch;
    coupling.mismatch.so.duration        = so_duration_mismatch;  
    
    coupling.so_total_sum                = so_all_amount;
    coupling.spindle_total_sum           = cell2mat(sum_coupling) + cell2mat(sum_mismatches);
    coupling.spindle_total_sum           = num2cell(coupling.spindle_total_sum);
      
    %%% Create percentage of match and mismatch
    empty_cells_index = cellfun('isempty', coupling.mismatch.sum);
    coupling.mismatch.sum(empty_cells_index) = {0};
    empty_cells_index = cellfun('isempty', coupling.match.sum);
    coupling.match.sum(empty_cells_index) = {0};

    spin_percentage_match              = (cell2mat(coupling.match.sum)./ cell2mat(coupling.spindle_total_sum))*100;
    coupling.match.spindle.percentage  = spin_percentage_match;
    so_percentage_match                = (cell2mat(coupling.match.sum)./ cell2mat(coupling.so_total_sum))*100;
    coupling.match.so.percentage       = so_percentage_match;
         
    %%% Save
    coupling_name = sprintf('coupling_sub%d.mat', subject_numbers(sub_nr));
    save(coupling_name, 'coupling')
    time_subj = toc;
    disp(['Finished subject ', num2str(sub_nr), ' in ', num2str(time_subj), ' seconds.' ])
           
   

end


