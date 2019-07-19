function [pca_all_mtx, comb_mtx_all, condition_idx_char, condition_idx_num] = ...
    pca_plot_last_drug(drug_test_sal_pow, drug_test_sal_pc, drug_post_sal_pow, drug_post_sal_pc, drug_pre_sal_pow, drug_pre_sal_pc, ...
    drug_test_MK_pow, drug_test_MK_pc, drug_post_MK_pow, drug_post_MK_pc, drug_pre_MK_pow, drug_pre_MK_pc,...
    drug_test_SC_pow, drug_test_SC_pc, drug_post_SC_pow, drug_post_SC_pc, drug_pre_SC_pow, drug_pre_SC_pc, bands, ret_nov_comb_mtx, condition_mtx_num)
%turn kevin's excel output matrices into a pca plot
%
%3 4
drop_subjs = [3 4];
keep_subjs = 1:size(drug_pre_sal_pow,1); keep_subjs = keep_subjs(~ismember(keep_subjs,drop_subjs));
num_sbjs = length(keep_subjs);

comb_mtx_sal = combine_mtxs(drug_test_sal_pow, drug_test_sal_pc, drug_post_sal_pow, drug_post_sal_pc, drug_pre_sal_pow, drug_pre_sal_pc, bands, keep_subjs);
comb_mtx_MK = combine_mtxs(drug_test_MK_pow, drug_test_MK_pc, drug_post_MK_pow, drug_post_MK_pc, drug_pre_MK_pow, drug_pre_MK_pc, bands, keep_subjs);
comb_mtx_SC = combine_mtxs(drug_test_SC_pow, drug_test_SC_pc, drug_post_SC_pow, drug_post_SC_pc, drug_pre_SC_pow, drug_pre_SC_pc, bands, keep_subjs);
comb_mtx_all = [comb_mtx_sal; comb_mtx_MK; comb_mtx_SC];
comb_mtx_all = [comb_mtx_all; ret_nov_comb_mtx]; 

%PREPROCESS TEST
%comb_mtx_all = comb_mtx_all./sum(comb_mtx_all,2);
%comb_mtx_all = log(comb_mtx_all);

%exclude conditions if desired
included_conditions = [6 7 9 10 12 13];

%included_conditions
%comb_mtx_all(ismember(condition_mtx_num, included_conditions),:) ;

%zscore
comb_mtx_all_z = zscore_mtx(comb_mtx_all);



%pca
if exist('comb_mtx_all_z', 'var')
    comb_mtx_all_pca = pca(comb_mtx_all_z');
else
    comb_mtx_all_pca = pca(comb_mtx_all');
end

pca_all_mtx = comb_mtx_all_pca;

%plot
figure; hold on

retnov_conditions{1,:} = 'retnov homecage ret';
retnov_conditions{2,:} = 'retnov homecage nov';
retnov_conditions{3,:} = 'retnov retrieval';
retnov_conditions{4,:} = 'retnov novel';

sizedrugmtx = size([comb_mtx_sal; comb_mtx_MK; comb_mtx_SC])
sizeretnovmtx = size(ret_nov_comb_mtx)
sizemtx = size(comb_mtx_all_pca)

num_drug_conditions = 9;
%plot ret_nov items, if they exist
%if size(comb_mtx_all_pca,1) > num_sbjs*num_drug_conditions
    ret_nov_rngs = [1 16; 17 48; 49 64; 65 96];
    ret_nov_rngs = ret_nov_rngs + repmat(num_sbjs*num_drug_conditions, size(ret_nov_rngs))
    colors2 = [...
        180 180 255;... %homecage (ret) (dark pastel blue)
        180 180 255;... %homecage (nov) (dark pastel blue)
        220 240 255;... %retrieval (light pastel blue)
        255 220 220;... %novelty (pastel pink)
        []]./255;
    for irn = 1 : size(ret_nov_rngs,1)
        
        condition_idx_char(ret_nov_rngs(irn,1): ret_nov_rngs(irn,2), :) = retnov_conditions(irn,:);
        condition_idx_num(ret_nov_rngs(irn,1): ret_nov_rngs(irn,2), :) = irn;
        handles_retnov{irn} = figfunction(comb_mtx_all_pca, ret_nov_rngs(irn,1), ret_nov_rngs(irn,2), colors2(irn,:), 0);
    end
%end

%plot drug items with black cirlces
drug_conditions{1,:} = 'drug sal test';
drug_conditions{2,:} = 'drug sal post';
drug_conditions{3,:} = 'drug sal pre';
drug_conditions{4,:} = 'drug MK test';
drug_conditions{5,:} = 'drug MK post';
drug_conditions{6,:} = 'drug MK pre';
drug_conditions{7,:} = 'drug SC test';
drug_conditions{8,:} = 'drug SC post';
drug_conditions{9,:} = 'drug SC pre';


colors = [0 200 220;... %sal test (light blue)
    0 120 220;... %sal post (mid blue)
    0 0 180;... %sal pre (dark blue)
    110 205 20;... %MK test (light green)
    100 160 50;... %MK post (mid green)
    0 128 0;... %MK pre (dark green)
    255 75 20;... %SC test (orange)
    210 35 10;... %SC post (mid red)
    130 20 5]; %SC pre (dark red)
colors = colors./255;
for icond = 1:num_drug_conditions

    condition_idx_char(1+num_sbjs*(icond-1) : num_sbjs*(icond), :) = drug_conditions(icond,:);
    condition_idx_num(1+num_sbjs*(icond-1) : num_sbjs*(icond), :) = irn+icond;
    handles_drug{icond} = figfunction(comb_mtx_all_pca, 1+num_sbjs*(icond-1), num_sbjs*(icond), colors(icond,:), 1);
    
end



handles_hold= [];

for ih = 1:length(handles_drug)
    handles_hold = [handles_hold handles_drug{ih}];
end
for ih = 1:length(handles_retnov)
    handles_hold = [handles_hold handles_retnov{ih}];
end

%legend(handles_hold, 'post SAL', 'post MK', 'post SC', 'test SAL', 'test MK', 'test SC')
legend(handles_hold, 'test SAL', 'post SAL', 'pre SAL', 'test MK', 'post MK', 'pre MK', 'test SC', 'post SC', 'pre SC', 'homecage (nov)', 'homecage (ret)', 'novelty', 'retrieval')

%'pre MK', 'pre SC','post SC',  'test SC')
%legend(handles_hold, 'pre SAL','post SAL','test SAL', 'pre MK','post MK','test MK')



axis square
xlim([min(comb_mtx_all_pca(:,1))*1.2 max(comb_mtx_all_pca(:,1))*1.2]); xlabel pc1
ylim([min(comb_mtx_all_pca(:,2))*1.2 max(comb_mtx_all_pca(:,2))*1.2]); ylabel pc2
zlim([min(comb_mtx_all_pca(:,3))*1.2 max(comb_mtx_all_pca(:,3))*1.2]); zlabel pc3
set(gca,'TickLength',[0, 0]);



%DISTANCES
%

%control preprocessing for distances
mtx_for_dists = comb_norm_mtx_all_z;

%preallocate ish
dists_out = [];

%combine homecage
%condition_idx(condition_idx==2) = 1;

%calculate iterative euclid distances
for ics= 1:size(mtx_for_dists,1)
    
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
    dists_out = [dists_out pdist([csesh; means_sans_csesh])'];
    
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


hc_and_ret = [ybar(3, ~isnan(ybar(3,:))) ybar(11, ~isnan(ybar(11,:)))];
hc_and_nov = [ybar(9, ~isnan(ybar(9,:))) ybar(17, ~isnan(ybar(17,:)))];
ret_and_nov = [ybar(14, ~isnan(ybar(14,:))) ybar(18, ~isnan(ybar(18,:)))];

%hc_and_ret = [ybar(11, ~isnan(ybar(11,:)))];
%hc_and_nov = [ybar(17, ~isnan(ybar(17,:)))];

figure; hold on;
%bar([mean(hc_and_ret) mean(hc_and_nov) mean(ret_and_nov)]) 
%errorbar([mean(hc_and_ret) mean(hc_and_nov) mean(ret_and_nov)], [std(hc_and_ret)/sqrt(length(hc_and_ret)) std(hc_and_nov)/sqrt(length(hc_and_nov)) std(ret_and_nov)/sqrt(length(ret_and_nov))], 'k.')
bar_out = [mean(hc_and_ret) mean(hc_and_nov) mean(ret_and_nov)]
bar([mean(hc_and_ret) mean(hc_and_nov)]) 
std_out_corrected = [std(hc_and_ret)/sqrt(length(hc_and_ret)) std(hc_and_nov)/sqrt(length(hc_and_nov)) std(ret_and_nov)/sqrt(length(ret_and_nov))]
errorbar([mean(hc_and_ret) mean(hc_and_nov)], [std(hc_and_ret)/sqrt(length(hc_and_ret)) std(hc_and_nov)/sqrt(length(hc_and_nov))], 'k.')

%[a b c d] = ttest2(hc_and_ret, hc_and_nov)

figure; hold on;  bar(mean(ret_and_nov)); errorbar(mean(ret_and_nov), std(ret_and_nov)/sqrt(length(ret_and_nov)))

[a b c d] = ttest(ret_and_nov, 0)
%}

%within functions
    function cmtx = rshp_comb(mtx1, mtx2, varargin)
        if nargin == 3
           bandz = varargin{1}; 
           keep_subjz = 1:size(mtx1,1);
        elseif nargin == 4
            bandz = varargin{1};
            keep_subjz = varargin{2};
        end
        cmtx = [mtx1(keep_subjz,:,bandz) mtx2(keep_subjz,:,bandz)];
        cmtx = reshape(cmtx, size(cmtx,1), size(cmtx,2)*size(cmtx,3));
        cmtx = cmtx(~isnan(sum(cmtx,2)),:);
    end

%combine in one matrix
    function cm = combine_mtxs(ip1, ip2, ip3, ip4, ip5, ip6, bands, varargin)
        
        if nargin == 8
            keep_subjz = varargin{1};
        else
            keep_subjz = 1:size(ip1,1);
        end
        
        c12 = rshp_comb(ip1, ip2, bands, keep_subjz);
        c34 = rshp_comb(ip3, ip4, bands, keep_subjz);
        c56 = rshp_comb(ip5, ip6, bands, keep_subjz);
        %c12 = (c12 - c56) ./ (c12 + c56); %normalize by pre
        %c34 = (c34 - c56) ./ (c34 + c56); %normalize by pre
        cm = [c12; c34; c56];
    end

%plot range of mtx rows
    function fighandle = figfunction(mtx, lo, high, color, encircled)
        
        dotsizes = [40 80];
        circlesizes = dotsizes./3.5;
        
        if encircled == 1
            plot3(mtx(lo:high, 1),...
            mtx(lo:high, 2),...
            mtx(lo:high, 3), '.', 'markersize', dotsizes(1), 'color', color);
        
            %small circles
            plot3(mtx(lo:high, 1),...
            mtx(lo:high, 2),...
            mtx(lo:high, 3), 'ko', 'markersize', circlesizes(1));
        
         fighandle = plot3(mean(mtx(lo:high, 1)),...
            mean(mtx(lo:high, 2)),...
            mean(mtx(lo:high, 3)), '.', 'markersize', dotsizes(2), 'color', color);
        
            %large circles
            plot3(mean(mtx(lo:high, 1)),...
            mean(mtx(lo:high, 2)),...
            mean(mtx(lo:high, 3)), 'ko', 'markersize', circlesizes(2));
        else
            plot3(mtx(lo:high, 1),...
            mtx(lo:high, 2),...
            mtx(lo:high, 3), '.', 'markersize', dotsizes(1), 'color', color);
         fighandle = plot3(mean(mtx(lo:high, 1)),...
            mean(mtx(lo:high, 2)),...
            mean(mtx(lo:high, 3)), '.', 'markersize', dotsizes(2), 'color', color);
        end
    end

end