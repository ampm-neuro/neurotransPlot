function byrat_catchall = plot_lfp_across_bands(ctx1, ctx2)


mean_catchall = nan(6, 3, 2);
std_catchall = nan(6, 3, 2);
byrat_catchall = nan(3, size(ctx1, 1), 6);
mean_mean = nan(1,6);

for band = 1:6
    
    [mc1, sc1, by_rat1] = lfp_by_freezing(ctx1, randperm(size(ctx1, 1)), band); ylim([0 2]); close
    [mc2, sc2, by_rat2] = lfp_by_freezing(ctx2, randperm(size(ctx2, 1)), band); ylim([0 2]); close
    
    mean_catchall(band, :, 1) = mc1; 
    mean_catchall(band, :, 2) = mc2;
    std_catchall(band, :, 1) = sc1; 
    std_catchall(band,:, 2) = sc2;
    byrat_catchall(:, :, band) = by_rat2 ./ (by_rat1+by_rat2);
    
    mean_mean(band) = nanmean(nanmean(nanmean(mean_catchall(band,:,:))));
    
end

figure
hold on
colors = [.1 .1 .1; .3 .3 .3; .6 .6 .6]; 
for pair = 1:3
    errorbar(mean_catchall(:,pair,1), std_catchall(:,pair,1)./sqrt(size(ctx1, 1)), 'Color', colors(pair,:),'linestyle','-')
    errorbar(mean_catchall(:,pair,2), std_catchall(:,pair,2)./sqrt(size(ctx1, 1)), 'Color', colors(pair,:),'linestyle', '--' )
end
xlim([.6 6.4])

%{
figure
hold on
colors = [.1 .1 .1; .3 .3 .3; .6 .6 .6]; 
for pair = 1:3

    errorbar(mean_catchall(:,pair,1)-mean_mean', std_catchall(:,pair,1)./sqrt(size(ctx1, 1)), 'Color', colors(pair,:),'linestyle','-')
    errorbar(mean_catchall(:,pair,2)-mean_mean', std_catchall(:,pair,2)./sqrt(size(ctx1, 1)), 'Color', colors(pair,:),'linestyle', '--' )
end
%}


%plot byrat figure
figure
hold on
ylim([-1 1])
byrat_mean = nanmean(byrat_catchall,2); byrat_mean = squeeze(byrat_mean);
byrat_std = nanstd(byrat_catchall, [], 2);  byrat_std = squeeze(byrat_std);
for pair = 1:3
    errorbar(byrat_mean(pair,:), byrat_std(pair,:)./sqrt(size(ctx1, 1)))
end
xlim([.6 6.4])
hold on; plot([.6 6.4], [.5 .5], 'k--')
ylim([.3 .69])

