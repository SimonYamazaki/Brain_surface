
fname = "/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/tutorial_data_20190918_1558/buckner_data/tutorial_subjs/004/surf/rh.pial";

[vertices, faces] = freesurfer_read_surf(fname);

%%
fn = "/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/GCF.csv";
T = readtable(fn).Var1;
colors_GC = T;

colors_GC(colors_GC<0) = 0;
colors_GC(colors_GC>0) = 1;

%% 
fn1 = "/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/MCF.csv";
T1 = readtable(fn1).Var1;
colors_MC = T1;

colors_MC(colors_MC<0) = 0;
colors_MC(colors_MC>0) = 1;

%% Gaussian curvature
%angle1
figure;
h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3))%,'FaceColor','green')
h.FaceVertexCData = colors_GC;
colormap spring
blueAndGreenColormap = [repmat([1,0,1], [1,1]) ; repmat([1,1,0], [1,1])]
colormap(blueAndGreenColormap) 

c = colorbar;
c.Location = 'north';
set(c,'Position',[0.15 0.92 0.75 0.03])

c.TicksMode = 'manual';
c.Ticks = [0,1];
c.TickLabels = {'Negative','Positive'};
c.FontSize = 12;


zoom(1.1)
view([-80,10,10])
set(gcf,'position',[0, 0, 500, 300]); 
saveas(gcf,'GC_bin_angle1.png')


%other angle
figure;
h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3))%,'FaceColor','green')
h.FaceVertexCData = colors_GC;
colormap spring
blueAndGreenColormap = [repmat([1,0,1], [1,1]) ; repmat([1,1,0], [1,1])]
colormap(blueAndGreenColormap) 

c = colorbar;
c.Location = 'north';
set(c,'Position',[0.15 0.92 0.75 0.03])

c.TicksMode = 'manual';
c.Ticks = [0,1];
c.TickLabels = {'Negative','Positive'};
c.FontSize = 12;

zoom(1.1)
view([50,10,10])
set(gcf,'position',[0, 0, 500, 300]); 
saveas(gcf,'GC_bin_angle2.png')


%zoom
figure;
h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3))%,'FaceColor','green')
h.FaceVertexCData = colors_GC;
colormap spring
blueAndGreenColormap = [repmat([1,0,1], [1,1]) ; repmat([1,1,0], [1,1])]
colormap(blueAndGreenColormap) 

c = colorbar;
c.Location = 'north';
set(c,'Position',[0.15 0.92 0.75 0.03])
set(c,'AxisLocation','in')

c.YTick = c.Limits;
c.YTickLabel = {'Negative','Positive'};
c.FontSize = 11;
c.FontWeight='bold';
c.Color = 'white';

ax = gca;  
ax.Clipping = 'off'; 

zoom(5)
view([50,0,230])
set(gcf,'position',[0, 0, 500, 300]); 
saveas(gcf,'GC_bin_zoom.png')

%% Mean Curvature

figure;
h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3))%,'FaceColor','green')
h.FaceVertexCData = colors_MC;
colormap spring
blueAndGreenColormap = [repmat([1,0,1], [1,1]) ; repmat([1,1,0], [1,1])]
colormap(blueAndGreenColormap) 

c = colorbar;
c.Location = 'north';
set(c,'Position',[0.15 0.92 0.75 0.03])

c.TicksMode = 'manual';
c.Ticks = [0,1];
c.TickLabels = {'Negative','Positive'};
c.FontSize = 10;

zoom(1.1)
view([-80,10,10])
set(gcf,'position',[0, 0, 500, 300]); 
saveas(gcf,'MC_bin_angle1.png')


%other angle 
figure;
h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3))%,'FaceColor','green')
h.FaceVertexCData = colors_MC;
colormap spring
blueAndGreenColormap = [repmat([1,0,1], [1,1]) ; repmat([1,1,0], [1,1])]
colormap(blueAndGreenColormap) 

c = colorbar;
c.Location = 'north';
set(c,'Position',[0.15 0.92 0.75 0.03])

c.TicksMode = 'manual';
c.Ticks = [0,1];
c.TickLabels = {'Negative','Positive'};
c.FontSize = 10;

zoom(1.1)
view([50,10,10])
set(gcf,'position',[0, 0, 500, 300]); 
saveas(gcf,'MC_bin_angle2.png')


%zoom
figure;
h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3))%,'FaceColor','green')
h.FaceVertexCData = colors_MC;
colormap spring
blueAndGreenColormap = [repmat([1,0,1], [1,1]) ; repmat([1,1,0], [1,1])]
colormap(blueAndGreenColormap) 

c = colorbar;
c.Location = 'north';
set(c,'Position',[0.15 0.92 0.75 0.03])
set(c,'AxisLocation','in')

c.YTick = c.Limits;
c.YTickLabel = {'Negative','Positive'};
c.FontSize = 11;
c.FontWeight='bold';
c.Color = 'white';

ax = gca;  
ax.Clipping = 'off'; 

zoom(5)
view([50,0,230])
set(gcf,'position',[0, 0, 500, 300]); 
saveas(gcf,'MC_bin_zoom.png')

%%
A=colormap('spring');