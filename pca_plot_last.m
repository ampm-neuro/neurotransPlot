function [hc_and_hc, hc_and_ret, hc_and_nov, ret_and_nov] = pca_plot_last(ret_pow, ret_pc, nov_pow, nov_pc, varargin)
%turn kevin's excel output matrices into a pca plot
%
% load('novel_retrieval.mat')
% try: pca_plot_last(ret_retrieval_power, ret_retrieval_PC, nov_novelty_power, nov_novelty_PC)

bands = 1:6;
%homecage inputs
if nargin > 4
    ret_pow_hc = varargin{1};
    ret_pc_hc = varargin{2};
    nov_pow_hc = varargin{3};
    nov_pc_hc = varargin{4};
    if length(varargin) == 5
        bands = varargin{5};
    end
end

%combine in one matrix
comb_norm_mtx_ret = rshp_comb(ret_pow, ret_pc, bands);

size(comb_norm_mtx_ret)

comb_norm_mtx_nov = rshp_comb(nov_pow, nov_pc, bands);

size(comb_norm_mtx_nov)

comb_norm_mtx_all = [comb_norm_mtx_ret; comb_norm_mtx_nov];

%add hc to top if exists
if exist('ret_pow_hc', 'var')
    comb_norm_mtx_ret_hc = rshp_comb(ret_pow_hc, ret_pc_hc, bands);
    comb_norm_mtx_nov_hc = rshp_comb(nov_pow_hc, nov_pc_hc, bands);
    comb_norm_mtx_all_hc = [comb_norm_mtx_ret_hc; comb_norm_mtx_nov_hc];
    comb_norm_mtx_all = [comb_norm_mtx_all_hc; comb_norm_mtx_all];
end


%log
%comb_norm_mtx_all = log(comb_norm_mtx_all);

%zscore
comb_norm_mtx_all_z = zscore_mtx(comb_norm_mtx_all);

%pca
if exist('comb_norm_mtx_all_z', 'var')
    comb_norm_mtx_all_pca = pca(comb_norm_mtx_all_z');
else
    comb_norm_mtx_all_pca = pca(comb_norm_mtx_all');
end

%plot
colors = [0 0 0;...
    0 50 150;...
    150 20 0;...
    80 80 80;...
    0 50 175;...
    210 30 0;...
    160 160 160;...
    0 50 220;...
    250 50 0]./255;

condition_idx = zeros(size(comb_norm_mtx_ret_hc(:,1)));
figure; hold on
if exist('comb_norm_mtx_all_hc', 'var')
    
    low = 1;
    hi = size(comb_norm_mtx_ret_hc,1);
    condition_idx(low:hi) = 1;
    plot3(comb_norm_mtx_all_pca(low:hi, 1), ...
        comb_norm_mtx_all_pca(low:hi, 2),...
        comb_norm_mtx_all_pca(low:hi, 3), '.', 'markersize', 40, 'color', colors(2,:));
    h1 = plot3(mean(comb_norm_mtx_all_pca(low:hi, 1)), ...
        mean(comb_norm_mtx_all_pca(low:hi, 2)),...
        mean(comb_norm_mtx_all_pca(low:hi, 3)), '.', 'markersize', 80, 'color', colors(2,:));
    
    
    low = hi+1;
    hi = hi + size(comb_norm_mtx_nov_hc,1);
    condition_idx(low:hi) = 2;
    plot3(comb_norm_mtx_all_pca(low:hi, 1),...
        comb_norm_mtx_all_pca(low:hi, 2),...
        comb_norm_mtx_all_pca(low:hi, 3), '.', 'markersize', 40, 'color', colors(3,:));
    h2 = plot3(mean(comb_norm_mtx_all_pca(low:hi, 1)),...
        mean(comb_norm_mtx_all_pca(low:hi, 2)),...
        mean(comb_norm_mtx_all_pca(low:hi, 3)), '.', 'markersize', 80, 'color', colors(3,:));
    
    
    low = hi+1;
    hi = hi + size(comb_norm_mtx_ret,1);
    condition_idx(low:hi) = 3;
    plot3(comb_norm_mtx_all_pca(low:hi, 1), ...
        comb_norm_mtx_all_pca(low:hi, 2),...
        comb_norm_mtx_all_pca(low:hi, 3), '.', 'markersize', 40, 'color', colors(8,:));
    h3 = plot3(mean(comb_norm_mtx_all_pca(low:hi, 1)), ...
        mean(comb_norm_mtx_all_pca(low:hi, 2)),...
        mean(comb_norm_mtx_all_pca(low:hi, 3)), '.', 'markersize', 80, 'color', colors(8,:));
    
    
    low = hi+1;
    hi = hi + size(comb_norm_mtx_nov,1);
    condition_idx(low:hi) = 4;
    plot3(comb_norm_mtx_all_pca(low:hi, 1),...
        comb_norm_mtx_all_pca(low:hi, 2),...
        comb_norm_mtx_all_pca(low:hi, 3), '.', 'markersize', 40, 'color', colors(6,:));
    h4 = plot3(mean(comb_norm_mtx_all_pca(low:hi, 1)),...
        mean(comb_norm_mtx_all_pca(low:hi, 2)),...
        mean(comb_norm_mtx_all_pca(low:hi, 3)), '.', 'markersize', 80, 'color', colors(6,:));

    legend([h1 h2 h3 h4], 'homecage (ret)', 'homecage (nov)', 'retrieval', 'novel')
else
    condition_idx(1:size(comb_norm_mtx_ret)) = 1;
    plot3(comb_norm_mtx_all_pca(1:size(comb_norm_mtx_ret,1),1), ...
        comb_norm_mtx_all_pca(1:size(comb_norm_mtx_ret,1),2),...
        comb_norm_mtx_all_pca(1:size(comb_norm_mtx_ret,1),3), 'b.', 'markersize', 30);
    condition_idx((size(comb_norm_mtx_ret,1)+1) : size(comb_norm_mtx_all_pca,1)) = 2;
    plot3(comb_norm_mtx_all_pca( (size(comb_norm_mtx_ret,1)+1) : end, 1),...
        comb_norm_mtx_all_pca( (size(comb_norm_mtx_ret,1)+1) : end, 2),...
        comb_norm_mtx_all_pca( (size(comb_norm_mtx_ret,1)+1) : end, 3), 'r.', 'markersize', 30);
end

xlabel x
ylabel y
zlabel z



%DISTANCES
%

%control preprocessing for distances
mtx_for_dists = comb_norm_mtx_all_z;

%preallocate ish
dists_out = [];


%condition_idx(condition_idx==2) = 1;

%calculate iterative euclid distances
for ics= 1:length(mtx_for_dists)
    
    %find current session
    csesh = mtx_for_dists(ics,:);
    
    %calculate means without the current sesh
    if exist('ret_pow_hc', 'var')
        means_sans_csesh = nan(4, size(mtx_for_dists,2));
        for im = 1:4
            type_idx = condition_idx == im;
            sesh_idx = ones(size(condition_idx)); sesh_idx(ics) = 1; sesh_idx = logical(sesh_idx);
            means_sans_csesh(im, :) = mean(mtx_for_dists(type_idx & sesh_idx,:));
        end
    else
        means_sans_csesh = nan(2, size(mtx_for_dists,2));
        for im = 1:2
            type_idx = condition_idx == im;
            sesh_idx = ones(size(condition_idx)); sesh_idx(ics) = 1; sesh_idx = logical(sesh_idx);
            means_sans_csesh(im, :) = mean(mtx_for_dists(type_idx & sesh_idx,:));
        end
    end
    
    %calculate and load distances
    dists_out = [dists_out (pdist([csesh; means_sans_csesh])./sqrt(length(csesh)))'];
    
end
    
%dist figure
%figure; hold on

ybar = nan(30, size(dists_out,2));
%std_prep = nan(30, 1);
%sqrt_num_samps = nan(1,30);
%distance from this (current) stage
for ib1 = 1:im %from above
    
    %distance to this stage
    for ib2 = 1:im
        xbar = 5*(ib1-1) + ib2;
        ybar(xbar, condition_idx == ib1) = (dists_out(ib2, condition_idx == ib1) - dists_out(ib1, condition_idx == ib1))./(dists_out(ib2, condition_idx == ib1) + dists_out(ib1, condition_idx == ib1));
        ybar_plot = mean(ybar(xbar, condition_idx == ib1));
        std_prep = std(ybar(xbar, condition_idx == ib1));
        
        %bar(xbar, ybar_plot)
        %errorbar(xbar, ybar_plot, std_prep/sqrt(sum(condition_idx == ib1)),'k.')
    end
    
end

%ybar is a weird matrix of arbitrary size(:,1) and size(:,2) of samples
%some rows are filled in above, along with columns that correspond to the
%conditions in condition_idx
%
%the first four rows are homecage_ret
%the next cluster of filled in rows are homecage_nov, then ret, then nov

hc_and_hc = [ybar(2, ~isnan(ybar(2,:))) ybar(6, ~isnan(ybar(6,:)))];
hc_and_ret = [ybar(3, ~isnan(ybar(3,:))) ybar(11, ~isnan(ybar(11,:)))];
hc_and_nov = [ybar(9, ~isnan(ybar(9,:))) ybar(17, ~isnan(ybar(17,:)))];
ret_and_nov = [ybar(14, ~isnan(ybar(14,:))) ybar(18, ~isnan(ybar(18,:)))];

%hc_and_ret = [ybar(11, ~isnan(ybar(11,:)))];
%hc_and_nov = [ybar(17, ~isnan(ybar(17,:)))];

figure; hold on;
%bar([mean(hc_and_ret) mean(hc_and_nov) mean(ret_and_nov)]) 
%errorbar([mean(hc_and_ret) mean(hc_and_nov) mean(ret_and_nov)], [std(hc_and_ret)/sqrt(length(hc_and_ret)) std(hc_and_nov)/sqrt(length(hc_and_nov)) std(ret_and_nov)/sqrt(length(ret_and_nov))], 'k.')
%bar_out = [mean(hc_and_ret) mean(hc_and_nov) mean(ret_and_nov)];
bar([mean(hc_and_hc) mean(hc_and_ret) mean(hc_and_nov) mean(ret_and_nov)]) 
%std_out_corrected = [std(hc_and_ret)/sqrt(length(hc_and_ret)) std(hc_and_nov)/sqrt(length(hc_and_nov)) std(ret_and_nov)/sqrt(length(ret_and_nov))];
errorbar([mean(hc_and_hc) mean(hc_and_ret) mean(hc_and_nov) mean(ret_and_nov)], [std(hc_and_hc)/sqrt(length(hc_and_hc)) std(hc_and_ret)/sqrt(length(hc_and_ret)) std(hc_and_nov)/sqrt(length(hc_and_nov)) std(ret_and_nov)/sqrt(length(ret_and_nov))], 'k.')

%[a b c d] = ttest2(hc_and_ret, hc_and_nov)

figure; hold on;  bar(mean(hc_and_hc)); errorbar(mean(hc_and_hc), std(hc_and_hc)/sqrt(length(hc_and_hc)))


[a b c d] = ttest(hc_and_hc, 0);


%within functions
    function cmtx = rshp_comb(mtx1, mtx2, varargin)
        if nargin == 3
           bandz = varargin{1}; 
        end
        cmtx = [mtx1(:,:,bandz) mtx2(:,:,bandz)];
        cmtx = reshape(cmtx, size(cmtx,1), size(cmtx,2)*size(cmtx,3));
        cmtx = cmtx(~isnan(sum(cmtx,2)),:);
    end
end