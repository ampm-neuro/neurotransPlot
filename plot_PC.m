function plot_PC(ctx1, ctx2)

%means, stds
mean_ctx1 = reshape(nanmean(ctx1), size(ctx1,2), size(ctx1,3))';
std_ctx1 = reshape(nanstd(ctx1), size(ctx1,2), size(ctx1,3))';
mean_ctx2 = reshape(nanmean(ctx2), size(ctx2,2), size(ctx2,3))';
std_ctx2 = reshape(nanstd(ctx2), size(ctx2,2), size(ctx2,3))';

%RSC-DH, RSC-ACC, RSC-DH
for col = 1:3
    figure; hold on; 

    errorbar(mean_ctx1(:,col), std_ctx1(:,col)./sqrt(squeeze(sum(~isnan(ctx1(:,col, :))))))
    errorbar(mean_ctx2(:,col), std_ctx2(:,col)./sqrt(squeeze(sum(~isnan(ctx2(:,col, :))))))
    
    title(num2str(col))
    ylim([0 2])
    xlim([.75 size(mean_ctx1,1)+.2])
end