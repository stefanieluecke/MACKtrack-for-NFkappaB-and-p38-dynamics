function h = colormapStack_multiple_NFkB(HMinput1,HMinput2,HMinput3,HMinput4,HMinput5,HMinput6)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Make stacked-colormap plot of cells
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



% Make figure (if not specified)
%{
if nargin < 7
    fig_handle = figure(gcf);% Create new figure, set properties
    set(fig_handle,'Color',[1 1 1],'InvertHardCopy', 'off','PaperPositionMode','auto')
end
%}

% Set figure/axes handles, 
%{
handles.figure1 = fig_handle;
handles.axes1 = axes('Parent', fig_handle);
handles.CellData = CellData;
handles.Measurement1 = measure1;
handles.Measurement2 = measure2;
handles.Options1 = options1;
handles.Options2 = options2;
%}

% Modify colormap so that time values outside of cell lifetime display as gray
mod_colormap = divergingmap(0:1/1023:1,[14 28 77]/255,[158 27 2]/255);%opens utilities function to make diverging color map with these specified rgb codes and the sepcified step size
mod_colormap(1,:) = [0.1 0.1 0.1]; %?

%attach measure1/CellData/Options to handles

%% Heatmap for ID1
% Make plot, setting axes and other options
%NFkB plot
subplot(1,6,1)

measure1DOWN = [zeros(1,size(HMinput1.var_nfkb_sorted,2));HMinput1.var_nfkb_sorted(1:end-1,:)];
measure1UP = [HMinput1.var_nfkb_sorted(2:end,:);zeros(1,size(HMinput1.var_nfkb_sorted,2))];
handles.h0 = imagesc(HMinput1.opt_nfkb.TimeBounds, [1 size(HMinput1.var_nfkb_sorted,1)],nan(size(HMinput1.var_nfkb_sorted)),HMinput1.opt_nfkb.MeasurementBounds);
hold on
handles.h1 = imagesc(HMinput1.opt_nfkb.TimeBounds, [1 size(HMinput1.var_nfkb_sorted,1)],measure1UP,HMinput1.opt_nfkb.MeasurementBounds);
handles.h2 = imagesc(HMinput1.opt_nfkb.TimeBounds, [1 size(HMinput1.var_nfkb_sorted,1)],measure1DOWN,HMinput1.opt_nfkb.MeasurementBounds);
handles.h3 = imagesc(HMinput1.opt_nfkb.TimeBounds, [1 size(HMinput1.var_nfkb_sorted,1)],HMinput1.var_nfkb_sorted,HMinput1.opt_nfkb.MeasurementBounds); 
hold off

colormap(mod_colormap)
c = colorbar('YTick',HMinput1.opt_nfkb.MeasurementTicks,'YTickLabel',HMinput1.opt_nfkb.MeasurementTickLabels);
set(c,'TickLength',0.003*ones(size(get(c,'TickLength'))))
ylabel(c,[HMinput1.opt_nfkb.Name],'FontSize',14);
xlabel('Time (h)','FontSize',14);
yticks([]);
xticks(HMinput1.opt_nfkb.TimeTicks);
title({HMinput1.title}, 'Interpreter', 'none')


%% Heatmap for ID2
% Make plot, setting axes and other options
%NFkB plot
subplot(1,6,2)

measure1DOWN = [zeros(1,size(HMinput2.var_nfkb_sorted,2));HMinput2.var_nfkb_sorted(1:end-1,:)];
measure1UP = [HMinput2.var_nfkb_sorted(2:end,:);zeros(1,size(HMinput2.var_nfkb_sorted,2))];
handles.h0 = imagesc(HMinput2.opt_nfkb.TimeBounds, [1 size(HMinput2.var_nfkb_sorted,1)],nan(size(HMinput2.var_nfkb_sorted)),HMinput2.opt_nfkb.MeasurementBounds);
hold on
handles.h1 = imagesc(HMinput2.opt_nfkb.TimeBounds, [1 size(HMinput2.var_nfkb_sorted,1)],measure1UP,HMinput2.opt_nfkb.MeasurementBounds);
handles.h2 = imagesc(HMinput2.opt_nfkb.TimeBounds, [1 size(HMinput2.var_nfkb_sorted,1)],measure1DOWN,HMinput2.opt_nfkb.MeasurementBounds);
handles.h3 = imagesc(HMinput2.opt_nfkb.TimeBounds, [1 size(HMinput2.var_nfkb_sorted,1)],HMinput2.var_nfkb_sorted,HMinput2.opt_nfkb.MeasurementBounds); 
hold off

colormap(mod_colormap)
c = colorbar('YTick',HMinput2.opt_nfkb.MeasurementTicks,'YTickLabel',HMinput2.opt_nfkb.MeasurementTickLabels);
set(c,'TickLength',0.003*ones(size(get(c,'TickLength'))))
ylabel(c,[HMinput2.opt_nfkb.Name],'FontSize',14);
xlabel('Time (h)','FontSize',14);
yticks([]);
xticks(HMinput2.opt_nfkb.TimeTicks);
title({HMinput2.title}, 'Interpreter', 'none')

%% Heatmap for ID3
% Make plot, setting axes and other options
%NFkB plot
subplot(1,6,3)

measure1DOWN = [zeros(1,size(HMinput3.var_nfkb_sorted,2));HMinput3.var_nfkb_sorted(1:end-1,:)];
measure1UP = [HMinput3.var_nfkb_sorted(2:end,:);zeros(1,size(HMinput3.var_nfkb_sorted,2))];
handles.h0 = imagesc(HMinput3.opt_nfkb.TimeBounds, [1 size(HMinput3.var_nfkb_sorted,1)],nan(size(HMinput3.var_nfkb_sorted)),HMinput3.opt_nfkb.MeasurementBounds);
hold on
handles.h1 = imagesc(HMinput3.opt_nfkb.TimeBounds, [1 size(HMinput3.var_nfkb_sorted,1)],measure1UP,HMinput3.opt_nfkb.MeasurementBounds);
handles.h2 = imagesc(HMinput3.opt_nfkb.TimeBounds, [1 size(HMinput3.var_nfkb_sorted,1)],measure1DOWN,HMinput3.opt_nfkb.MeasurementBounds);
handles.h3 = imagesc(HMinput3.opt_nfkb.TimeBounds, [1 size(HMinput3.var_nfkb_sorted,1)],HMinput3.var_nfkb_sorted,HMinput3.opt_nfkb.MeasurementBounds); 
hold off

colormap(mod_colormap)
c = colorbar('YTick',HMinput3.opt_nfkb.MeasurementTicks,'YTickLabel',HMinput3.opt_nfkb.MeasurementTickLabels);
set(c,'TickLength',0.003*ones(size(get(c,'TickLength'))))
ylabel(c,[HMinput3.opt_nfkb.Name],'FontSize',14);
xlabel('Time (h)','FontSize',14);
yticks([]);
xticks(HMinput3.opt_nfkb.TimeTicks);
title({HMinput3.title}, 'Interpreter', 'none')


%% Heatmap for ID4
% Make plot, setting axes and other options
%NFkB plot
subplot(1,6,4)

measure1DOWN = [zeros(1,size(HMinput4.var_nfkb_sorted,2));HMinput4.var_nfkb_sorted(1:end-1,:)];
measure1UP = [HMinput4.var_nfkb_sorted(2:end,:);zeros(1,size(HMinput4.var_nfkb_sorted,2))];
handles.h0 = imagesc(HMinput4.opt_nfkb.TimeBounds, [1 size(HMinput4.var_nfkb_sorted,1)],nan(size(HMinput4.var_nfkb_sorted)),HMinput4.opt_nfkb.MeasurementBounds);
hold on
handles.h1 = imagesc(HMinput4.opt_nfkb.TimeBounds, [1 size(HMinput4.var_nfkb_sorted,1)],measure1UP,HMinput4.opt_nfkb.MeasurementBounds);
handles.h2 = imagesc(HMinput4.opt_nfkb.TimeBounds, [1 size(HMinput4.var_nfkb_sorted,1)],measure1DOWN,HMinput4.opt_nfkb.MeasurementBounds);
handles.h3 = imagesc(HMinput4.opt_nfkb.TimeBounds, [1 size(HMinput4.var_nfkb_sorted,1)],HMinput4.var_nfkb_sorted,HMinput4.opt_nfkb.MeasurementBounds); 
hold off

colormap(mod_colormap)
c = colorbar('YTick',HMinput4.opt_nfkb.MeasurementTicks,'YTickLabel',HMinput4.opt_nfkb.MeasurementTickLabels);
set(c,'TickLength',0.003*ones(size(get(c,'TickLength'))))
ylabel(c,[HMinput4.opt_nfkb.Name],'FontSize',14);
xlabel('Time (h)','FontSize',14);
yticks([]);
xticks(HMinput4.opt_nfkb.TimeTicks);
title({HMinput4.title}, 'Interpreter', 'none')

%% Heatmap for ID5
% Make plot, setting axes and other options
%NFkB plot
subplot(1,6,5)

measure1DOWN = [zeros(1,size(HMinput5.var_nfkb_sorted,2));HMinput5.var_nfkb_sorted(1:end-1,:)];
measure1UP = [HMinput5.var_nfkb_sorted(2:end,:);zeros(1,size(HMinput5.var_nfkb_sorted,2))];
handles.h0 = imagesc(HMinput5.opt_nfkb.TimeBounds, [1 size(HMinput5.var_nfkb_sorted,1)],nan(size(HMinput5.var_nfkb_sorted)),HMinput5.opt_nfkb.MeasurementBounds);
hold on
handles.h1 = imagesc(HMinput5.opt_nfkb.TimeBounds, [1 size(HMinput5.var_nfkb_sorted,1)],measure1UP,HMinput5.opt_nfkb.MeasurementBounds);
handles.h2 = imagesc(HMinput5.opt_nfkb.TimeBounds, [1 size(HMinput5.var_nfkb_sorted,1)],measure1DOWN,HMinput5.opt_nfkb.MeasurementBounds);
handles.h3 = imagesc(HMinput5.opt_nfkb.TimeBounds, [1 size(HMinput5.var_nfkb_sorted,1)],HMinput5.var_nfkb_sorted,HMinput5.opt_nfkb.MeasurementBounds); 
hold off

colormap(mod_colormap)
c = colorbar('YTick',HMinput5.opt_nfkb.MeasurementTicks,'YTickLabel',HMinput5.opt_nfkb.MeasurementTickLabels);
set(c,'TickLength',0.003*ones(size(get(c,'TickLength'))))
ylabel(c,[HMinput5.opt_nfkb.Name],'FontSize',14);
xlabel('Time (h)','FontSize',14);
yticks([]);
xticks(HMinput5.opt_nfkb.TimeTicks);
title({HMinput5.title}, 'Interpreter', 'none')

%% Heatmap for ID6
% Make plot, setting axes and other options
%NFkB plot
subplot(1,6,6)

measure1DOWN = [zeros(1,size(HMinput6.var_nfkb_sorted,2));HMinput6.var_nfkb_sorted(1:end-1,:)];
measure1UP = [HMinput6.var_nfkb_sorted(2:end,:);zeros(1,size(HMinput6.var_nfkb_sorted,2))];
handles.h0 = imagesc(HMinput6.opt_nfkb.TimeBounds, [1 size(HMinput6.var_nfkb_sorted,1)],nan(size(HMinput6.var_nfkb_sorted)),HMinput6.opt_nfkb.MeasurementBounds);
hold on
handles.h1 = imagesc(HMinput6.opt_nfkb.TimeBounds, [1 size(HMinput6.var_nfkb_sorted,1)],measure1UP,HMinput6.opt_nfkb.MeasurementBounds);
handles.h2 = imagesc(HMinput6.opt_nfkb.TimeBounds, [1 size(HMinput6.var_nfkb_sorted,1)],measure1DOWN,HMinput6.opt_nfkb.MeasurementBounds);
handles.h3 = imagesc(HMinput6.opt_nfkb.TimeBounds, [1 size(HMinput6.var_nfkb_sorted,1)],HMinput6.var_nfkb_sorted,HMinput6.opt_nfkb.MeasurementBounds); 
hold off

colormap(mod_colormap)
c = colorbar('YTick',HMinput6.opt_nfkb.MeasurementTicks,'YTickLabel',HMinput6.opt_nfkb.MeasurementTickLabels);
set(c,'TickLength',0.003*ones(size(get(c,'TickLength'))))
ylabel(c,[HMinput6.opt_nfkb.Name],'FontSize',14);
xlabel('Time (h)','FontSize',14);
yticks([]);
xticks(HMinput6.opt_nfkb.TimeTicks);
title({HMinput6.title}, 'Interpreter', 'none')

