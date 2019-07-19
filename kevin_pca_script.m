%file prompt
%
[filename, pathname] = uigetfile('*.*');

%load file
%
[mtx, ~, raw] = xlsread([pathname filename]);
[grouping_variables, ~, gv_idx] = unique(raw(2:end,2));

%report columns and rows (make sure these are right!)
%
number_of_mice = size(mtx, 1)/length(grouping_variables)
grouping_variables = grouping_variables
number_of_DVs = size(mtx, 2)

%Remove rats with missing data (using mean to detect NaNs)
%
idx_hold = reshape(1:size(mtx,1), sum(gv_idx==1), 1,...
    length(unique(gv_idx)));
mtx_hold = mean(mean(permute(reshape(mtx',...
    [size(mtx,2), size(mtx,1)/length(unique(gv_idx)),...
    length(unique(gv_idx))]), [2,1,3]), 3), 2);
idx_hold = idx_hold(~isnan(mtx_hold),:,:);
if length(idx_hold) < size(mtx,1) %correct and report
    mice_removed = unique(mice(setdiff(1:length(mice), idx_hold)))
    mtx = mtx(idx_hold(:), :);
    gv_idx = gv_idx(idx_hold(:));
end


%zscore data matrix
%
mtx = mtx - repmat(nanmean(mtx), size(mtx,1), 1);
mtx_z = mtx./nanstd(mtx);
    %remove columns without variance
    mtx_z = mtx_z(:,~isnan(mtx_z(1,:)) & ~isinf(mtx_z(1,:)));

%calculate pca
%
[pca_mtx, ~, ~, ~, pct_explained_var] = pca(mtx_z');

%select to-be-plotted principal components
%
PCA_cols_to_plot = [1 2 3];
plot_vars = cell(length(PCA_cols_to_plot), 1);
handles = cell(length(grouping_variables), 1);
if any(~ismember(PCA_cols_to_plot, 1:size(pca_mtx,2)))
    PCA_cols_to_plot = PCA_cols_to_plot
    size_pca_mtx = size(pca_mtx)
    error(['You are requesting more dimensions'...
        ' than you have principal components'])
end

%check dimensions and plot
%
figure;hold on;
ml_colors = get(gca,'ColorOrder');
sds = 40; %small dot size
lds = 120; %large dot size (means)

%iterate through grouping variables
for igv = unique(gv_idx)'

    %get plotting variables
    for ipc = 1:length(PCA_cols_to_plot)
        plot_vars{ipc} = pca_mtx(gv_idx==igv, PCA_cols_to_plot(ipc));
    end

    %plot plotting variables
    if length(PCA_cols_to_plot) == 1 %if only 1 pc
        
        handles{igv} = plot(plot_vars{1}, '-', 'linewidth', 2,...
            'color', ml_colors(igv, :));
        
    elseif length(PCA_cols_to_plot) == 2 %if 2 pcs
        
        handles{igv} = plot(plot_vars{1}, plot_vars{2},...
            '.', 'Markersize', sds, 'color', ml_colors(igv, :)); %data points
        plot(mean(plot_vars{1}), mean(plot_vars{2}),...
            '.', 'Markersize', lds, 'color', ml_colors(igv, :)) %means
        
    elseif length(PCA_cols_to_plot) == 3 %if 3 pcs
        
        handles{igv} = plot3(plot_vars{1}, plot_vars{2}, plot_vars{3}, '.',...
            'Markersize', sds, 'color', ml_colors(igv, :)); %data points
        plot3(mean(plot_vars{1}), mean(plot_vars{2}),...
            mean(plot_vars{3}), '.', 'Markersize', lds, 'color',...
            ml_colors(igv, :)) %means
        
    else % if other number of pcs
        error(['check PCA_cols_to_plot specification. '...
            'must be a vector of only 1 2 or 3 numbers.'])
    end
end

%figure beautification
%
%legend
handles_hold= [];
for ih = 1:length(handles)
    handles_hold = [handles_hold handles{ih}];
end
legend(handles_hold, grouping_variables,'Location','NorthEastOutside')

%axis labels
xlabel(['Principal Component ' num2str(PCA_cols_to_plot(1))])
if length(PCA_cols_to_plot) == 2 %if 2 pcs
    ylabel(['Principal Component ' num2str(PCA_cols_to_plot(2))])
elseif length(PCA_cols_to_plot) == 3 %if 3 pcs
    ylabel(['Principal Component ' num2str(PCA_cols_to_plot(2))])
    zlabel(['Principal Component ' num2str(PCA_cols_to_plot(3))])
end

%misc
set(gca,'TickLength',[0, 0]); 
box off
axis square

