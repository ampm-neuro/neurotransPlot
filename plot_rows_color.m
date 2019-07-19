function plot_rows_color(pcs, dims, row_chunk, varargin)
%plot matrix with different colors for row chunks

if nargin == 4
    colors = varargin{1};
else
    figure
    colors = get(gca,'colororder'); 
    colors = [colors;colors;colors];
    close
end


figure; hold on

clust_cell = cell((size(pcs,1)/row_chunk),1);

idx_min = 1;
for i = 1:(size(pcs,1)/row_chunk)
    
    idx_max = row_chunk*i;
    
    
    idx_rng = idx_min:idx_max;
    %idx_rng = idx_rng(good_rows(idx_rng));
    
    
    clust_cell{i} = pcs(idx_rng, dims);
    
    if numel(dims) == 2
        plot(pcs(idx_rng, dims(1)), pcs(idx_rng, dims(2)), '.', 'MarkerSize', 50, 'Color', colors(i,:))
        plot(nanmean(pcs(idx_rng, dims(1))), nanmean(pcs(idx_rng, dims(2))), '.', 'MarkerSize', 150, 'Color', colors(i,:))
    elseif numel(dims) == 3
        plot3(pcs(idx_rng, dims(1)), pcs(idx_rng, dims(2)), pcs(idx_rng, dims(3)), '.', 'MarkerSize', 50, 'Color', colors(i,:))
        plot3(nanmean(pcs(idx_rng, dims(1))), nanmean(pcs(idx_rng, dims(2))), nanmean(pcs(idx_rng, dims(3))), '.', 'MarkerSize', 150, 'Color', colors(i,:))
    else
        error('incorrect dims input')
        
    end
    
    idx_min = idx_max+1;



end





xlabel(['pc' num2str(dims(1))])
ylabel(['pc' num2str(dims(2))])
zlabel(['pc' num2str(dims(3))])
