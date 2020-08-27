function out_mean = get_mean(in_vals, in_chan, in_chan_vec)
% example: spin_duration = get_mean(spinso_coupling.test_duration_seconds, channel_info_subj(1,:));
    out_mean = zeros(1,length(in_chan_vec));
    for k = 1:size(in_chan_vec,2)
        mask = strcmp(in_chan, in_chan_vec(k));
        val = mean(in_vals(mask));
        out_mean(k) = val;
    end
end   
