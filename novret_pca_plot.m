figure; hold on

colors = [0,0.547,0.8910;...
    0.89,0.30,0.098;...
    0.929,0.694,0.125;...
    0.494,0.184,0.556;...
    0.466,0.674,0.188;...
    0.301,0.745,0.933;...
    0.635,0.078,0.184];

dims = [3 4 5];

rng = cell(4,1);
rng{1} = 1:33;
rng{2} = 34:50;
rng{3} = 51:83;
rng{4} = 84:100;

plot3(novret_plot_mtx(rng{1},dims(1)), novret_plot_mtx(rng{1},dims(2)), novret_plot_mtx(rng{1},dims(3)), '.', 'markersize', 50, 'color', colors(1,:)./1.45)
plot3(novret_plot_mtx(rng{2},dims(1)), novret_plot_mtx(rng{2},dims(2)), novret_plot_mtx(rng{2},dims(3)), '.', 'markersize', 50, 'color', colors(2,:)./1.45)
plot3(novret_plot_mtx(rng{3},dims(1)), novret_plot_mtx(rng{3},dims(2)), novret_plot_mtx(rng{3},dims(3)), '.', 'markersize', 50, 'color', colors(1,:))
plot3(novret_plot_mtx(rng{4},dims(1)), novret_plot_mtx(rng{4},dims(2)), novret_plot_mtx(rng{4},dims(3)), '.', 'markersize', 50, 'color', colors(2,:))

plot3(nanmean(novret_plot_mtx(rng{1},dims(1))), nanmean(novret_plot_mtx(rng{1},dims(2))), nanmean(novret_plot_mtx(rng{1},dims(3))), '.', 'markersize', 120, 'color', colors(1,:)./1.45);
plot3(nanmean(novret_plot_mtx(rng{2},dims(1))), nanmean(novret_plot_mtx(rng{2},dims(2))), nanmean(novret_plot_mtx(rng{2},dims(3))), '.', 'markersize', 120, 'color', colors(2,:)./1.45);
plot3(nanmean(novret_plot_mtx(rng{3},dims(1))), nanmean(novret_plot_mtx(rng{3},dims(2))), nanmean(novret_plot_mtx(rng{3},dims(3))), '.', 'markersize', 120, 'color', colors(1,:));
plot3(nanmean(novret_plot_mtx(rng{4},dims(1))), nanmean(novret_plot_mtx(rng{4},dims(2))), nanmean(novret_plot_mtx(rng{4},dims(3))), '.', 'markersize', 120, 'color', colors(2,:));

box off; set(gca,'TickLength',[0, 0]);

xlabel('pc1')
ylabel pc2
zlabel('pc3')

if isequal(dims, [1 2 3])
xlim([.075 .13])
ylim([-.22 .34])
zlim([-.27 .3])
end

%{
hold_ = novret_plot_mtx(84:100, dims);
freez_rng = freezing./max(freezing); freez_rng(isnan(freez_rng)) = .5;
for ih = 1:size(hold_,1)
    plot3(hold_(ih,dims(1)), hold_(ih,dims(2)), hold_(ih,dims(3)), '.', 'markersize', 50, 'color', colors(2,:).*freez_rng(ih))
end
%}