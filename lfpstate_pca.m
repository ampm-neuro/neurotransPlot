function [dist_mtx, plot_mtx] = lfpstate_pca(power_cell, pc_cell, bands)
%plots pca of lfp band data
%
%input cells contain multiple 3d matrices, each from one recording
%condition. cells correspond between inputs. E.g., power_cell(1) comes from
%the same recording session as pc_cell(1).
%
%bands is a vector containing any 1:6, specifying which lfp bands to
%include in the computation

%items per condition (number of subjects)
ipc = size(power_cell{1},1);

%modify each rat's (individually) power and then coherence
%pull out the rats slices from each condition
%{
for irat = 1:ipc
    
    %pull out power slices
    power_hold_comb = [];
    for iphc = 1:length(power_cell)
        power_catch = power_cell{iphc}(irat,:,:);
        power_hold_comb = [power_hold_comb; power_catch(:)];
    end
    %zscore power slices
    power_hold_comb = power_hold_comb-repmat(nanmean(power_hold_comb),size(power_hold_comb));
    power_hold_comb = power_hold_comb./repmat(nanstd(power_hold_comb),size(power_hold_comb)); 
    %put back power slices
    power_hold_comb = reshape(power_hold_comb, size(power_hold_comb,1)/length(power_cell), length(power_cell)); %two columns
    
    for icol = 1:size(power_hold_comb,2)
        power_cell{icol}(irat,:,:) = reshape(power_hold_comb(:,icol), size(power_cell{iphc}(irat,:,:)));
    end
end
%}


%merge cells into single matrix (multiply coherence and power items)
%
powerPC_mtx = [];
for ic = 1:length(power_cell)

    %simple combine
    combine_pwrpc = [pc_cell{ic}(:,:,bands) power_cell{ic}(:,:,bands)];
    combine_pwrpc = reshape(combine_pwrpc, size(combine_pwrpc,1), size(combine_pwrpc,2)*size(combine_pwrpc,3));%reshape
    
    %load
    powerPC_mtx = [powerPC_mtx; combine_pwrpc];

end


figure; imagesc(powerPC_mtx)

%rows sans nans
rsn = sum(isnan(powerPC_mtx),2)==0;
num_removed_rows = sum(~rsn)

%remove rows with nans
powerPC_mtx_sansnans = powerPC_mtx(rsn, :);

figure; imagesc(powerPC_mtx_sansnans); title start

%log
%powerPC_mtx = powerPC_mtx;
%powerPC_mtx = log(powerPC_mtx);

%zscore
powerPC_mtx_sansnans = zscore_mtx(powerPC_mtx_sansnans);
figure; imagesc(powerPC_mtx_sansnans); title zscore


dist_mtx = nan(size(powerPC_mtx,1), size(powerPC_mtx_sansnans,2)); 
dist_mtx(rsn, :) = powerPC_mtx_sansnans;


%{
cluster_a = dist_mtx(1:ipc, :);
cluster_b = dist_mtx(ipc+1:2*ipc, :);

%dists
zc_dist = zscore_2cluster_dist(cluster_a, cluster_b)
try
    m_dist = mahal_2cluster_dist(cluster_a, cluster_b)
catch
end
%}

%pca
pca_mtx = pca(powerPC_mtx_sansnans');


%replace nan rows for plotting
plot_mtx = nan(size(powerPC_mtx,1), size(pca_mtx,2)); 
plot_mtx(rsn, :) = pca_mtx;



%colors
colors = [0 0 0;...
    0 50 150;...
    150 20 0;...
    80 80 80;...
    0 50 175;...
    210 30 0;...
    160 160 160;...
    0 50 220;...
    250 50 0]./255;


%plot
plot_rows_color(plot_mtx, [1 2 3], ipc, colors)


end