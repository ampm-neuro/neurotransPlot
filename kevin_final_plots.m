function [pw_dists_out, xlabel_cell, comb_mtx_all_pca, condition_mtx_char] = kevin_final_plots(comb_all_mtx, condition_mtx_num, condition_mtx_char)
%cuts up all lfp info and makes plots

%exclude undesired subjects
drug_outliers = [];
drug_outliers = [3 4];
exclude = [];
for id = unique(drug_outliers)
    exclude = [id:7:7*9  exclude];
end
exclude = unique(exclude);
comb_all_mtx = comb_all_mtx(setdiff(1:size(comb_all_mtx,1), exclude),:);
condition_mtx_num = condition_mtx_num(setdiff(1:size(condition_mtx_num,1), exclude),:);
condition_mtx_char = condition_mtx_char(setdiff(1:size(condition_mtx_char,1), exclude),:);


%exclude undesired conditions
%included_conditions = [1 2 3 4]; %RetNov
%included_conditions = [6 7 9 10 12 13]; %PrePost
%included_conditions = [5 6 8 9 11 12]; %PostTest
included_conditions = unique(condition_mtx_num); %all
included_idx = ismember(condition_mtx_num, included_conditions);
comb_all_mtx = comb_all_mtx(included_idx,:);
condition_mtx_char = condition_mtx_char(included_idx); [~, unq_idx] = unique(condition_mtx_char, 'first'); included_conditions = condition_mtx_char(sort(unq_idx));
condition_mtx_num = condition_mtx_num(included_idx);

%zscore
comb_mtx_all_z = zscore_mtx(comb_all_mtx);

%pca
comb_mtx_all_pca = pca(comb_mtx_all_z');

%plot pca
if length(included_conditions) > 6
    colors = zeros(length(included_conditions), 3);
elseif length(included_conditions) == 6
    colors = colorfun([1 3 4 6 7 9]);
elseif length(included_conditions) == 4
    colors = colorfun([3 6 1 4]);
else
    error('colors incorrect')
end

fighandle = cell(1,length(unique(condition_mtx_num)));
figure; hold on
count = 0;
for ic = unique(condition_mtx_num)'
    count = count+1;
    idx = condition_mtx_num == ic;        
    fighandle{count} = figfunction(comb_mtx_all_pca, idx, colors(count,:), 1);
    %fighandle{count} = figfunction(comb_mtx_all_pca, idx, colors(1,:), 1); %1 color
end

handles_hold = [];
for ih = 1:length(fighandle)
    handles_hold = [handles_hold fighandle{ih}];
end

legend(handles_hold, included_conditions)
axis square
xlim([min(comb_mtx_all_pca(:,1))*1.2 max(comb_mtx_all_pca(:,1))*1.2]); xlabel pc1
ylim([min(comb_mtx_all_pca(:,2))*1.2 max(comb_mtx_all_pca(:,2))*1.2]); ylabel pc2
zlim([min(comb_mtx_all_pca(:,3))*1.2 max(comb_mtx_all_pca(:,3))*1.2]); zlabel pc3
set(gca,'TickLength',[0, 0]);

%dists
pw_dists_out = [];
pw_dists_out_lengths = [];
xlabel_cell = [];
comparisons = [1000 1000];
counter = 0;
for ic1 = unique(condition_mtx_num)'
    for ic2 = unique(condition_mtx_num)'
        
        %skip identity
        if isequal(ic1,ic2) || max(sum(ismember(comparisons, [ic1 ic2]),2))==2
            continue
        else
            counter = counter +1;
        end
        
        pw_dists_out(:,counter) = nan(size(comb_all_mtx,1), 1);

        dists_out_hold = twoclustdist(comb_all_mtx(condition_mtx_num == ic1, :), comb_all_mtx(condition_mtx_num == ic2, :));
        pw_dists_out_lengths(counter) = length(dists_out_hold);
        pw_dists_out(1:pw_dists_out_lengths(counter), counter) = dists_out_hold;
        
        comparisons(counter+1,:) = [ic1 ic2];
    
ahold = strjoin([included_conditions(unique(condition_mtx_num)==ic1) ' ' included_conditions(unique(condition_mtx_num)==ic2)]);
xlabel_cell{counter} = ahold;
        
    end
end


%plotdists
figure; hold on
bar(nanmean(pw_dists_out));
set(gca,'TickLength',[0, 0]);
errorbar(nanmean(pw_dists_out),nanstd(pw_dists_out)./sqrt(pw_dists_out_lengths),'k.')
set(gca,'FontSize',4)
xticks(1:length(nanmean(pw_dists_out)))
xticklabels(xlabel_cell)



%iterative euclidian dists
    function dists_out = twoclustdist(mtx1, mtx2)
        
        dists_out = nan(size([mtx1;mtx2],1),1);
        
        %iterate through first matrix
        for ics1= 1:size(mtx1,1)
            
            %current
            cs1 = mtx1(ics1,:);
            
            %means sans current
            mtx1_mean = mean(mtx1(setdiff(1:size(mtx1,1), ics1),:));
            mtx2_mean = mean(mtx2);
            
            %dists
            dists_hold = pdist([cs1; mtx1_mean; mtx2_mean])./sqrt(size([mtx1;mtx2],2));
            dists_out(ics1) = (dists_hold(2) - dists_hold(1)) / (dists_hold(2) + dists_hold(1));
            
        end

        %iterate through second matrix
        for ics2 = 1:size(mtx2,1)
            
            %current
            cs2 = mtx2(ics2,:);
            
            %means sans current
            mtx1_mean = mean(mtx1);
            mtx2_mean = mean(mtx2(setdiff(1:size(mtx2,1), ics2),:));
            
            %dists
            dists_hold = pdist([cs2; mtx2_mean; mtx1_mean])./sqrt(size([mtx1;mtx2],2));
            dists_out(ics1+ics2) = (dists_hold(2) - dists_hold(1)) / (dists_hold(2) + dists_hold(1));
            
        end

    end

%colors
    function colors = colorfun(colorsel)
        
        colors = [0 200 220;... %light blue
                0 120 220;... %mid blue
                0 0 180;... %dark blue
                110 215 20;... %light green
                100 160 50;... %mid green
                0 128 0;... %dark green
                255 55 20;... %orange
                210 35 10;... %mid red
                170 20 5]; %dark red
            
        colors = colors(colorsel,:)./255; 
    end

%plot range of mtx rows
    function fighandle = figfunction(mtx, idx, color, encircled)
        
        dotsizes = [40 80];
        circlesizes = dotsizes./3.5;
        
        if encircled == 1
            plot3(mtx(idx, 1),...
            mtx(idx, 2),...
            mtx(idx, 3), '.', 'markersize', dotsizes(1), 'color', color);
        
            %small circles
            plot3(mtx(idx, 1),...
            mtx(idx, 2),...
            mtx(idx, 3), 'ko', 'markersize', circlesizes(1));
        
         fighandle = plot3(mean(mtx(idx, 1)),...
            mean(mtx(idx, 2)),...
            mean(mtx(idx, 3)), '.', 'markersize', dotsizes(2), 'color', color);
        
            %large circles
            plot3(mean(mtx(idx, 1)),...
            mean(mtx(idx, 2)),...
            mean(mtx(idx, 3)), 'ko', 'markersize', circlesizes(2));
        else
            plot3(mtx(idx, 1),...
            mtx(idx, 2),...
            mtx(idx, 3), '.', 'markersize', dotsizes(1), 'color', color);
         fighandle = plot3(mean(mtx(idx, 1)),...
            mean(mtx(idx, 2)),...
            mean(mtx(idx, 3)), '.', 'markersize', dotsizes(2), 'color', color);
        end
    end

end