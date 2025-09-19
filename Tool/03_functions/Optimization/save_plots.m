figs = findall(0, 'Type', 'figure');
fmin_fig = findobj(figs, 'Name', 'Optimization Plot Function');
fPSO_fig = findobj(figs, 'Name', 'particleswarm');


fig_opt=figure('Position',[3841 -318 1536 747],'Color',[1,1,1]);
set(fig_opt,'defaultAxesTickLabelInterpreter','latex');  
set(fig_opt,'defaulttextinterpreter','latex');
set(fig_opt,'defaultLegendInterpreter','latex');
sgtitle(['Optimization Test: ',T_OPT_ii.test_name{:}])

ax_=subplot(1,2,1);
hold on
grid on


title(ax_,['Particle Swarm search (PSO) - Best Value: ', num2str(min(findobj(fPSO_fig.Children,'Type','Axes').Children.YData))])

set(ax_,'LineWidth',1.5,...
    'Fontsize',14);
xlabel(ax_,'Iteration')
ylabel(ax_,'Best Modify $H_{2}$ norm value')

plot(ax_,findobj(fPSO_fig.Children,'Type','Axes').Children.XData,findobj(fPSO_fig.Children,'Type','Axes').Children.YData,...
    'Color',"k")
scatter(ax_,findobj(fPSO_fig.Children,'Type','Axes').Children.XData,findobj(fPSO_fig.Children,'Type','Axes').Children.YData,...
    50,"k", 'filled')


ax_=subplot(1,2,2);
hold on
grid on
title(ax_,['Local search (fmin) - Best Value: ', num2str(min(findobj(fmin_fig.Children,'Type','Axes').Children.YData))])
set(ax_,'LineWidth',1.5,...
    'Fontsize',14);
xlabel(ax_,'Iteration')
ylabel(ax_,'Best Modify $H_{2}$ norm value')

plot(ax_,findobj(fmin_fig.Children,'Type','Axes').Children.XData,findobj(fmin_fig.Children,'Type','Axes').Children.YData,...
     'Color',"k")
scatter(ax_,findobj(fmin_fig.Children,'Type','Axes').Children.XData,findobj(fmin_fig.Children,'Type','Axes').Children.YData,...
    50,"k", 'filled')

savefig(fig_opt,[folder_results,'Optimization_process.fig'])
exportgraphics(fig_opt,[folder_results,'Optimization_process.png'], 'Resolution', 300)


