function [means_catch, stds_catch, by_rat] = lfp_by_freezing(ret_lfp, frz_times, varargin)
%compares lfp during retrieval between subjects with different freezing
%times


n_grps = 1;%number of groups
n_subjs = sum(~isnan(frz_times));%number_subj
subj_bin_size = floor(n_subjs/n_grps);%min subjects per group

%freezing sort index
[frz_times_sort, frz_idx] = sort(frz_times);



%divide up ret_lfp by sorted subject groups
%e.g., ret_retrieval (rats, pairs, bands)
%
%preallocate groupings
cell_grps = cell(n_grps, 1); %low to high sort

%which pairs and bands do we want?
pairs = [1 2 3];

if nargin > 2
    bands = varargin{1};
else
    bands = [1 2 3 4 5 6];
    %bands = [2 3]; %theta
    %bands = [1 4 5 6];%~theta
    %bands = [1];
end

ct_low = 1;
for ig = 1:n_grps
    if ig < n_grps
        ct_hi = ct_low+subj_bin_size-1;
        frx_idx_local = frz_idx(ct_low:ct_hi); %grab relevant subjects
        ct_low = ct_hi+1;
    elseif ig == n_grps
         frx_idx_local = frz_idx(ct_low:end);
    end
    
    %pairs, bands, rats
    cell_grps{ig} = permute(ret_lfp(frx_idx_local, pairs, bands), [2, 3, 1]);%load cells
    
end

%compute desired averages and stds
means_catch = [];
stds_catch = [];
for ic = 1:length(cell_grps)
    
    %output individual rat values
    by_rat = squeeze(cell_grps{ic});
    
    %combine bands within rats
    mean_local = nanmean(cell_grps{ic},2);
    
    %calculate mean and std across rats
    means_catch = [means_catch nanmean(mean_local,3)];
    stds_catch = [stds_catch nanstd(mean_local, [], 3)];
    
end

%plot
figure; hold on
bar(means_catch')
errorbar(means_catch', stds_catch'./sqrt(n_subjs), 'k.')



%scatterplot
%{
sp = permute(ret_lfp(:, pairs, bands), [2, 3, 1]);
sp_comb = nanmean(sp,2); sp_comb = permute(sp_comb, [3, 1, 2]);

[R, p] = fit_line(frz_times(~isnan(frz_times)),sp_comb(~isnan(frz_times),1))
title pairs1
[R, p] = fit_line(frz_times(~isnan(frz_times)),sp_comb(~isnan(frz_times),2))
title pairs2
%}





end