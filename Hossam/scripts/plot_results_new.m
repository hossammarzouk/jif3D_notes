%% load gravity responses
%%
projectPath = uigetdir;
plotRMS(projectPath)
drawnow
iter = inputdlg('Iterantion number',"write an iteration number");
iter = str2num(iter{1});

%%

projectPath="Y:\jointInversion\EW_2024\m3_crossGradient";
iter = 12;
%%
projectPath="D:\jif3d_JointInversion\EW_profile_Egypt\m7_joint_MI_7err";
iter = 45;

%%
projectPath = "R:\jointInversion\EW_Egypt\m7_increaseMTweight";
iter=47;
%%
plotRMS(projectPath)

%%
projectPath="Y:\jointInversion\EW_2024\m7_joint_MI_7err";
iter = 46;
%%
load("D:\PhD_Thesis\scripts\ch_06_inversion\sobh_moho_lab.mat")
[mtModel.X mtModel.Y ]= deg2utm(mtModel.lat,mtModel.long,'36 R');

EqclippedUTM = readtable('D:\Egypt_EW_all\GIS\Eq_clipped_UTM_usedME.csv')

disp("Loading..")
grav = 1;mt = 1; mag=0;

% load("DataInfo.mat")
cmap = load("D:\Software\utils\crameri colormap\CrameriColourMaps7.0.mat");
[~,model_name,~] = fileparts(projectPath);

mode = "UTM"; % "model" ; "UTM"
% mtFile.name   = fullfile(projectPath,['result' num2str(iter) '.mt.inv.nc']);
if mt
    mtFile.name   = dir(fullfile(projectPath,['*' num2str(iter) '.mt.inv.nc']));
    mtFile.name=fullfile(mtFile.name.folder,mtFile.name.name);
    mtFile.Northing = ncread(mtFile.name,'Northing');
    mtFile.Easting  = ncread(mtFile.name,'Easting');
    mtFile.Depth    = ncread(mtFile.name,'Depth');
    mtFile.origin   =   [ncread(mtFile.name,'Northing_Origin'), ncread(mtFile.name,'Easting_Origin'),...
        ncread(mtFile.name,'Depth_Origin')];
    mtFile.resistivity = log10(1/ncread(mtFile.name,'Conductivity'));
    [x y z] = meshgrid(mtFile.Easting-MT.yshift,mtFile.Northing-MT.xshift,mtFile.Depth);
    disp(['MT file loaded: ✅ ' mtFile.name ])
end
% if mag
%     magFile.name   = dir(fullfile(projectPath,['*' num2str(iter) '.mt.inv.nc']));
%     magFile.name=fullfile(magFile.name.folder,mtFile.name.name);
%     magFile.Northing = ncread(magFile.name,'Northing');
%     magFile.Easting  = ncread(magFile.name,'Easting');
%     magFile.Depth    = ncread(magFile.name,'Depth');
%     magFile.origin   =   [ncread(magFile.name,'Northing_Origin'), ncread(magFile.name,'Easting_Origin'),...
%         ncread(mtFile.name,'Depth_Origin')];
%     mtFile.resistivity = log10(1/ncread(mtFile.name,'Conductivity'));
%     [x y z] = meshgrid(mtFile.Easting-MT.yshift,mtFile.Northing-MT.xshift,magFile.Depth);
%     disp(['MT file loaded: ✅ ' mtFile.name ])
% end
if grav
    gravFile.name   = dir(fullfile(projectPath,['*t' num2str(iter) '.grav.inv.nc']));
    gravFile.name=fullfile(gravFile.name.folder,gravFile.name.name);
    gravFile.Northing = ncread(gravFile.name,'Northing');
    gravFile.Easting  = ncread(gravFile.name,'Easting');
    gravFile.Depth    = ncread(gravFile.name,'Depth');
    gravFile.origin   =   [ncread(gravFile.name,'Northing_Origin'), ncread(gravFile.name,'Easting_Origin'),...
        ncread(gravFile.name,'Depth_Origin')];
    gravFile.Density = ncread(gravFile.name,'Density');
    % gravFile.Density = ncread(gravFile.name,'Conductivity');

    [x y z] = meshgrid(gravFile.Easting-gravity.yshift,gravFile.Northing-gravity.xshift,gravFile.Depth);
    disp(['Gravity file loaded: ✅ ' gravFile.name])
end

% read respones file

if grav
    gravFileResp.name   = dir(fullfile(projectPath,'*.inv_sgd.nc'));
    gravFileResp.name=fullfile(gravFileResp.name.folder,gravFileResp.name.name);
    gravFileResp.Northing = ncread(gravFileResp.name,'MeasPosX');
    gravFileResp.Easting  = ncread(gravFileResp.name,'MeasPosY');
    gravFileResp.Depth    = ncread(gravFileResp.name,'MeasPosZ');
    gravFileResp.Density = ncread(gravFileResp.name,'Scalar_gravity');
    % [x y z] = meshgrid(gravFileResp.Easting-gravity.yshift,gravFileResp.Northing-gravity.xshift,gravFileResp.Depth);
    disp(['Gravity response file loaded: ✅ ' gravFileResp.name])
end
%% depth slice

ax={};
fig = figure('Position',[ 752         700        900 503]);
% fig = figure('Position',[ 752         700        900 503],'SizeChangedFcn',@fig_resize);

% fig = figure('Position',[ 752         700        900         503]);


% depthSlices = [0 1000 5000 10000 20000 30000 40000 50000 80000 100000 120000 150000];
depthSlices = [1000  5000 8000 12000  16000 20000 25000 35000 40000 50000];

tile = tiledlayout(numel(depthSlices),2);
tile.TileSpacing = "tight";
tile.Title.String = model_name;
counter=1;
for i = 1:numel(depthSlices)
    % disp(['Plotting... Model: ' num2str(i_model) ' Depth:' num2str(i) '/' num2str(numel(depthSlices))]);
    % if isequal(i_model,1)
    % ax{i} = nexttile(counter);
    % if isequal(i,1)
    %     title(model_name);
    % end
% 
% [xr yr] = rotate_points(EqclippedUTM.X,EqclippedUTM.Y,-MT.rotation,MT.originUTM(1),MT.originUTM(2));
% xr = xr+ MT.originUTM(1);
% yr = yr + MT.originUTM(2);

 % scatter( EqclippedUTM.X,EqclippedUTM.Depth_km_*1000,    EqclippedUTM.ML*10,'filled','MarkerFaceColor','b')

 % scatter( xr,EqclippedUTM.Depth_km_*1000,    EqclippedUTM.ML*10,'filled','MarkerFaceColor','r')

    % prepareAxis(ax{i,i_model},mode)

    if isequal(mode,"UTM")
        % axis 'auto xy'

        buffer=50000;
        % zlim([-2000 165000])
        % ylim([min(sx2{i_model})-buffer max(sx2{i_model})+buffer+buffer])
        % xlim([min(sy2{i_model})-buffer max(sy2{i_model})+buffer])

        % hold on
        if grav
            ax{i,1} = nexttile;
            sl = slice(x,y,z,gravFile.Density,[],[],depthSlices(i));
            % rotate(sl,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0])
            hold on
            % p1 = plot(MT.Yrot,MT.Xrot,'o','Color','k','MarkerFaceColor','k','MarkerSize',2);

            xlim([min(gravFile.Easting)-gravity.yshift max(gravFile.Easting)-gravity.yshift]);
            ylim([min(gravFile.Northing)-gravity.xshift max(gravFile.Northing)-gravity.xshift]);
            title(['Depth ' num2str(depthSlices(i)/1000) ' km'],'VerticalAlignment','baseline'   );
            colormap(ax{i,1},cmap.roma);colorbar
            caxis([-50 50]);
            view(2);
            shading flat

            plotshpFiles_UTM()

                set(gca,'SortMethod','childorder','XTickLabel',[])
        end

        if mt
            ax{i,2} = nexttile;
            sl = slice(x,y,z,mtFile.resistivity,[],[],depthSlices(i));
            % rotate(sl,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0]);hold on
            % p1 = plot(MT.Yrot,MT.Xrot,'o','Color','k','MarkerFaceColor','k','MarkerSize',2);
            % xlim([min(MT.Yrot)-buffer max(MT.Yrot)+buffer]);
            % ylim([min(MT.Xrot)-buffer max(MT.Xrot)+buffer]);
                    title(['Depth ' num2str(depthSlices(i)/1000) ' km'],'VerticalAlignment','baseline'   );
 colormap(ax{i,2},jet);colorbar
    caxis([0 4]);
    view(2);
    shading flat

    plotshpFiles_UTM()


                set(gca,'SortMethod','childorder','XTickLabel',[])

        end
        % rotate(sl,[0 0 1],rotation,[y0,x0,0])

    else
      if grav
            ax{i,1} = nexttile;
            sl = slice(x,y,z,gravFile.Density,[],[],depthSlices(i));
            % rotate(sl,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0])
            hold on
            % p1 = plot(MT.Yrot,MT.Xrot,'o','Color','k','MarkerFaceColor','k','MarkerSize',2);

            % xlim([min(gravFile.Easting)-gravity.yshift max(gravFile.Easting)-gravity.yshift]);
            % ylim([min(gravFile.Northing)-gravity.xshift max(gravFile.Northing)-gravity.xshift]);
            title(['Depth ' num2str(depthSlices(i)/1000) ' km'],'VerticalAlignment','baseline'   );
            colormap(ax{i,1},cmap.roma);colorbar
            caxis([-50 50]);
            view(2);
            shading flat

            % plotshpFiles_UTM()
                set(gca,'SortMethod','childorder','XTickLabel',[])
        end

        if mt
            ax{i,2} = nexttile;
            sl = slice(x,y,z,mtFile.resistivity,[],[],depthSlices(i));
            % rotate(sl,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0]);hold on
            % p1 = plot(MT.Yrot,MT.Xrot,'o','Color','k','MarkerFaceColor','k','MarkerSize',2);
            % xlim([min(MT.Yrot)-buffer max(MT.Yrot)+buffer]);
            % ylim([min(MT.Xrot)-buffer max(MT.Xrot)+buffer]);
                    title(['Depth ' num2str(depthSlices(i)/1000) ' km'],'VerticalAlignment','baseline'   );
 colormap(ax{i,2},jet);colorbar
    caxis([0 4]);
    view(2);
    shading flat

    % plotshpFiles_UTM()
                set(gca,'SortMethod','childorder','XTickLabel',[])

        end
    end
    % colormap(flip(pmkmp(50,'CubicL')))
    hold on
    % xlim([min(gravFile.Easting)-gravity.yshift max(gravFile.Easting)-gravity.yshift]);
    % ylim([min(gravFile.Northing)-gravity.xshift max(gravFile.Northing)-gravity.xshift]);

    % m1 =  plot3(sy2{i_model},sx2{i_model}, sz2{i_model}-1000,'marker','o','markersize',2,'color','k','linestyle','none','markerfacecolor','k');
    % rotate(m1,[0 0 1],rotation,[y0,x0,0])

    % colormap(flip(jet));
   

    % % plot hotspots
    % rt = text(ax{i,i_model}.XLim(2)-0.2,ax{i}.YLim(1)+.2,['Depth ' num2str(zb(idx)/1000,'%.0f') ' km'], ...
    %     'HorizontalAlignment','right', ...
    %     'VerticalAlignment','bottom', ...
    %     'FontSize',8,'Color','black');
    % rectangle(ax{i,i_model},'Position',[rt.Position(1)-1.5,rt.Position(2),2,.5],'Curvature',0.2,'FaceColor','white');
    % rt = text(ax{i,i_model}.XLim(2)-0.2,ax{i}.YLim(1)+.2,['Depth ' num2str(zb(idx)/1000,'%.0f') ' km'], ...
    %     'HorizontalAlignment','right', ...
    %     'VerticalAlignment','bottom', ...
    %     'FontSize',8,'Color','black');
    % set(gca,'SortMethod','childorder'),
    counter=counter+1;


end
% 
% cb = colorbar('LineWidth',1.5,...
%     'FontWeight','bold',...
%     'FontSize',12);
% cb.Label.Interpreter = "latex";
% cb.Label.String = 'Density [mGal]';
% % cb.Label.String = '$\log_{10}(Resistivity[\Omega.m])$';
% cb.Layout.Tile = 'south';

% Link axes po
grav_ax = [];  res_ax=[];
for iax = 1:size(ax,1)
    grav_ax = [grav_ax ax{iax,1}];
    res_ax = [res_ax ax{iax,2}];
   
    % if iseven(iax)
    %     set(ax{iax},"YTickLabel",[],"YLabel",[]);
    % end
    % if ~(ismember(iax, [size(ax,2), size(ax,2)-1]))
    %     set(ax{iax},"XTickLabel",[],"XLabel",[]);
    % end
end
%
linkaxes(res_ax, "xyz");
linkprop(res_ax,"CLim");
linkaxes(grav_ax, "xyz");
linkprop(grav_ax,"CLim");
%%
exportgraphics(gcf,fullfile(projectPath,'depthSlices_09_new.png'),'resolution',300)
%%
plotshpFiles_UTM()
%% plot resuidal and response diffrede
figure
tile = tiledlayout(1,4)
ax1 = nexttile;
scatter(gravity.Y+gravity.yshift,gravity.X+gravity.xshift,[],gravity.Data.*1e5,'filled');
colormap(cmap.roma);colorbar;
caxis([-50 50])
ax2 = nexttile;
scatter(gravFileResp.Easting+gravity.yshift,gravFileResp.Northing+gravity.xshift,[],gravFileResp.Density.*1e5,'filled');
colormap(cmap.roma);colorbar;
caxis([-50 50])

ax3 = nexttile;
scatter(gravFileResp.Easting+gravity.yshift,gravFileResp.Northing+gravity.xshift,[],(gravFileResp.Density-gravity.Data).*1e5,'filled');
colormap(ax3,jet);colorbar;
caxis([-50 50])

ax4 = nexttile;
histogram(gravFileResp.Density-gravity.Data.*1e5);
% colormap(ax4,jet);colorbar;
% caxis([-50 50])












%% Depth slice
sliceDepth = 1000;
fig = figure; fig.Position = [ 1000         563        1249         775];
tile = tiledlayout('flow');
ax1 = nexttile;
hold on
set(ax1,'FontWeight', 'bold','FontSize',10,'Box','on','LineWidth',1.5,'ZDir','reverse');
s1 = slice(ax1, x,y,z,mtFile.resistivity,[],[],sliceDepth) % resistivity
rotate(s1,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0])
p1 = plot(MT.Yrot,MT.Xrot,'v','Color','k')
colormap(ax1,flip(pmkmp(50,'CubicL')));
colorbar
shading flat
axis tight
ax1.Title.String = ['Resistivity model at ' num2str(sliceDepth) ' m'];
%
% ax2 = nexttile;
% hold on
% set(ax2,'FontWeight', 'bold','FontSize',10,'Box','on','LineWidth',1.5,'ZDir','reverse');
% slice(ax2, x,y,z,gravFile.Density,[],[],sliceDepth)
% colormap(ax2,"viridis");
% colorbar
% ax2.Title.String = ['Density model at ' num2str(sliceDepth) ' m'];

shading flat
axis tight
% Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
% ylim([2700000,3100000])
view(2)
%%
plotshpFiles_UTM()
%% slice only on stations locations
% hold on
roi = [];
% exclude station from profile
ch2 = string(strtrim(MT.sitesName));
ch_ind = find(strtrim(ch2)~='soh1');

% C = sortrows([MT.Yshifted(ch_ind)-MT.yshift MT.Xshifted(ch_ind)-MT.xshift ],'ascend');
C = sortrows([MT.Y(ch_ind) MT.X(ch_ind)],'ascend');

roi.Position  = [double(C(:,1)) double(C(:,2))];

% [s1 s2] = sort(sy2(ch_ind));
% roi.Position  = [sy2(s2) sx2(s2)];

% plot(roi.Position(:,1),roi.Position(:,2),'LineWidth',2,'Color','red');
% plot(sy2,sx2,'LineWidth',2,'LineStyle','none','Marker','o');

buffer = 10000; % buffer out of profile
tmp  = interp1([roi.Position(1,1),roi.Position(2,1)], [roi.Position(1,2),roi.Position(2,2)], ...
    roi.Position(1,1)-buffer,"linear","extrap");
roi.Position(1,:)    =  [roi.Position(1,1)-buffer, tmp];

tmp  = interp1([roi.Position(end-1,1),roi.Position(end,1)], [roi.Position(end-1,2),roi.Position(end,2)], ...
    roi.Position(end,1)+buffer,"linear","extrap");
roi.Position(end,:)    =  [roi.Position(end,1)+buffer, tmp];


% hold on
spacing = 0.1;
xq = .1 : spacing:numel(roi.Position(:,1));
vq = interp1(roi.Position(:,2),xq);

xq = .1 : spacing :numel(roi.Position(:,2));
vq2 = interp1(roi.Position(:,1),xq);

[xd,zd]=meshgrid(vq,squeeze(z(1,1,:)));
[yd,zd]=meshgrid(vq2,squeeze(z(1,1,:)));
% arbityry slice plot
fig = figure; fig.Position = [ 1000         563        1249         775];
tile = tiledlayout('flow');
ax = {};
if grav
ax{1} = nexttile;
hold on
target = ax{end}
set(target,'FontWeight', 'bold','FontSize',10,'Box','on','LineWidth',1.5,'ZDir','reverse','YDir','reverse');
slRes= slice(x,y,z,gravFile.Density,yd,xd,zd);

rotate(slRes,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0])
% colormap(target,flip(pmkmp(50,'CubicL')));
colormap(ax{1},cmap.roma);
colorbar
shading flat
% axis tight
target.Title.String = ['Resistivity model Joint inversion'];
view(0,0)
caxis([-50 50])

p1 = plot(MT.Yrot,MT.Xrot,'v','Color','k');
% rotate(p1,[0 0 1],15,[y0,x0,0])
plotshpFiles_UTM()
            scatter3( EqclippedUTM.X,EqclippedUTM.Y,EqclippedUTM.Depth_km_*1000,    EqclippedUTM.ML*10,'filled','MarkerFaceColor','r')

p = plot3(mtModel.X,mtModel.Y,mtModel.moho*1000,'k','LineWidth',2);
rotate(p,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0])

p = plot3(mtModel.X,mtModel.Y,mtModel.LAB*1000,'k','LineWidth',2);
rotate(p,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0])

  xlim([min(MT.Yrot)-buffer max(MT.Yrot)+buffer]);
            ylim([min(MT.Xrot)-buffer max(MT.Xrot)+buffer]);
end
if mt
ax{2} = nexttile;
hold on
target = ax{end}
set(target,'FontWeight', 'bold','FontSize',10,'Box','on','LineWidth',1.5,'ZDir','reverse','YDir','reverse');
slRes= slice(x,y,z,mtFile.resistivity,yd,xd,zd);
rotate(slRes,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0])
% colormap(target,flip(pmkmp(50,'CubicL')));
colormap(ax{2},flip(jet))
colorbar
shading flat
% axis tight
target.Title.String = ['Resistivity model Joint inversion'];
view(0,0)
p1 = plot(MT.Yrot,MT.Xrot,'v','Color','k');
% rotate(p1,[0 0 1],15,[y0,x0,0])
plotshpFiles_UTM()
            scatter3( EqclippedUTM.X,EqclippedUTM.Y,EqclippedUTM.Depth_km_*1000,    EqclippedUTM.ML*10,'filled','MarkerFaceColor','r')

p = plot3(mtModel.X,mtModel.Y,mtModel.moho*1000,'k','LineWidth',2);
rotate(p,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0])

p = plot3(mtModel.X,mtModel.Y,mtModel.LAB*1000,'k','LineWidth',2);
rotate(p,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0])

  xlim([min(MT.Yrot)-buffer max(MT.Yrot)+buffer]);
            ylim([min(MT.Xrot)-buffer max(MT.Xrot)+buffer]);
caxis([0 4])
% slRes= slice(x,y,z,gravFile.Density,yd,xd,zd);
                set(gca,'SortMethod','childorder')


% ax{2} = nexttile;
% plot3(target,sy2,sx2, sz,'marker','v','markersize',4,'color','k','linestyle','none','markerfacecolor','k');
% plot3(target,goldLoc(1),goldLoc(2), 0,'marker','diamond','markersize',4,'color','red','linestyle','none','markerfacecolor','k');
%
% target = ax{end}
% hold on
% set(target,'FontWeight', 'bold','FontSize',10,'Box','on','LineWidth',1.5,'ZDir','reverse','YDir','reverse');
% slDen= slice(target,x,y,z,gravFile.Density,yd,xd,zd);
% % colormap(target,"viridis");
% cmap = cbrewer2('Spectral');
% colormap(target,flip(cmap));
% colorbar
% target.Title.String = ['Density model Joint inversion' ];
% shading flat
% view(0,0)
% plot3(target,sy2,sx2, sz,'marker','v','markersize',4,'color','k','linestyle','none','markerfacecolor','k');
% plot3(target,goldLoc(1),goldLoc(2), 0,'marker','diamond','markersize',4,'color','red','linestyle','none','markerfacecolor','k');
%
end
Link = linkprop([ax{1}, ax{2}],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
%% plot on 2d surfaces
fig = figure; fig.Position = [ 1000         563        1249         775];
tile = tiledlayout(2,1);
ax = {};
ax{1} = nexttile;
hold on
target = ax{end}
set(target,'FontWeight', 'bold','FontSize',10,'Box','on','LineWidth',1.5,'ZDir','reverse','YDir','reverse');
t1 = interp3(x,y,z,gravFile.Density,yd,xd,zd);
slRes = surf(yd,zd,t1)
hold on
% [con1 con] = contour(yd,zd,10.^t1,[190 190 ],'ShowText','on','LineColor',"black",...
    % "LineWidth",1.5);
% [con1 con] = contour(yd,zd,10.^t1,[100:100:400],'ShowText','on','LineColor',"white",    "LineWidth",0.5);


% rotate(slRes,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0])
% colormap(target,flip(pmkmp(50,'CubicL')));
colormap(ax{1},cmap.roma);
colorbar
shading flat
% axis tight
target.Title.String = ['Density model Joint inversion'];
view(2)
caxis([-50 50])

plot(yd(:),mtModel.moho*1000,'w','LineWidth',2);
plot(yd(:),mtModel.LAB*1000,'w','LineWidth',2);

% p1 = plot(MT.Yrot,MT.Xrot,'v','Color','k');
% rotate(p1,[0 0 1],15,[y0,x0,0])
% plotshpFiles_UTM()
  xlim([min(MT.Yrot)-buffer max(MT.Yrot)+buffer]);
            % ylim([min(MT.Xrot)-buffer max(MT.Xrot)+buffer]);
ax{2} = nexttile;
hold on
target = ax{end}
set(target,'FontWeight', 'bold','FontSize',10,'Box','on','LineWidth',1.5,'ZDir','reverse','YDir','reverse');
t1 = interp3(x,y,z,mtFile.resistivity,yd,xd,zd);
slRes = surf(yd,zd,t1)
hold on
% rotate(slRes,[0 0 1],-MT.rotation,[MT.originUTM(1),MT.originUTM(2),0])
% colormap(target,flip(pmkmp(50,'CubicL')));
colormap(ax{2},flip(jet))
colorbar
shading flat
% axis tight
target.Title.String = ['Resistivity model Joint inversion'];
view(2)
% p1 = plot(MT.Yrot,MT.Xrot,'v','Color','k');
% rotate(p1,[0 0 1],15,[y0,x0,0])
% plotshpFiles_UTM()
  xlim([min(MT.Yrot)-buffer max(MT.Yrot)+buffer]);
            % ylim([min(MT.Xrot)-buffer max(MT.Xrot)+buffer]);
caxis([0 4])
% slRes= slice(x,y,z,gravFile.Density,yd,xd,zd);

plot(yd(:),mtModel.moho*1000,'k','LineWidth',2);
plot(yd(:),mtModel.LAB*1000,'k','LineWidth',2);

%% profile 2d for respnses
% d1 = interp2(gravity.Y,gravity.X,gravity.Data ,yd,xd);
% F1  = scatteredInterpolant(gravity.Y-yshift,gravity.X-xshift,gravity.Data );
% F2  = scatteredInterpolant(gravFileResp.Easting,gravFileResp.Northing,gravFileResp.Density );
figure;
F1  = scatteredInterpolant(gravity.Y,gravity.X,gravity.Data );
F2  = scatteredInterpolant(gravFileResp.Easting-gravity.yshift,gravFileResp.Northing-gravity.xshift,gravFileResp.Density );
d1 = F1(yd(:),xd(:));
d2 = F2(yd(:),xd(:));
% [xd,zd]=meshgrid(vq,squeeze(z(1,1,:)));
% [yd,zd]=meshgrid(vq2,squeeze(z(1,1,:)));
% plot(yd(:),detrend(d1*1e5,'omitnan'),'LineStyle','-',LineWidth=2);
% plot(yd(:),d1*1e5,'LineStyle','-',LineWidth=2);

hold on
% plot(yd(:),detrend(d2*1e5,'omitnan'),'LineStyle','--',LineWidth=2,Color='r')
plot(yd(:),detrend(d2*1e5,'omitnan'  ) ,'LineStyle','--',LineWidth=2,Color='r')

err = ones(numel(d1),1);
% errorbar(yd(:),d1*1e5,err)
legend({"Obs";"Calc"})
shadedErrorBar(yd(:),detrend(d1*1e5,"omitnan"),err,'lineProps','-b')
%% plot parameters realations
% figure
tile = tiledlayout(1,2);tile.TileSpacing="tight";
cm = cbrewer2('Spectral',10);
maxDepth = 100000;
minDepth = 1;

nexttile
z_indx = find(zd(:,1)>minDepth & zd(:,1)<maxDepth);
zd2 = zd(z_indx,:);
slgrav= interp3(x,y,z,gravFile.Density,yd(z_indx,:),xd(z_indx,:),zd(z_indx,:));
slRes= interp3(x,y,z,mtFile.resistivity,yd(z_indx,:),xd(z_indx,:),zd(z_indx,:));
scatter(slgrav(:),slRes(:),[],reshape(zd(z_indx,:),1,[]),'Filled');hold on
colormap(cm); caxis([minDepth maxDepth])
xlabel("Density");ylabel("Resistivity");title("profile location");box on

nexttile
d = permute(gravFile.Depth(z_indx),[3 2 1]);
d= repmat(d,[size(gravFile.Density,1),size(gravFile.Density,2),1]);
scatter(reshape(gravFile.Density(:,:,z_indx),1,[]),reshape(mtFile.resistivity(:,:,z_indx),1,[]),[],d(:),'Filled');hold on
colormap(cm); caxis([minDepth maxDepth])
xlabel("Density");title("All Model");box on

cb = colorbar('LineWidth',1.5,...
    'FontWeight','bold',...
    'FontSize',12);
% cb.Label.Interpreter = "latex";
cb.Label.String = 'Depth (km)';
% cb.Label.String = '$\log_{10}(Resistivity[\Omega.m])$';
cb.Layout.Tile = 'south';cb.TickLabelsMode="auto";
cb.TickLabels  = cellfun(@(x) str2num(x)*10 ,cb.TickLabels  ,'UniformOutput' ,false);
% leg = {};
% for i = 1:size(slRes,1)
% plot(slgrav(i,:),slRes(i,:),'+','MarkerFaceColor','auto');hold on
% leg{end+1} = num2str(mtFile.Depth(i));
% legend(leg);
% end

%% counts (2d hsitograms)
h = histogram2(slgrav(1:10,:),slRes(1:10,:),'FaceColor','flat','EdgeColor','none','NumBins',[200 200]);view(2)
caxis([1 0100])
%%
function plotshpFiles_UTM()
% Load shapefiles
hold on
river = shaperead('D:\Egypt_EW_all\GIS\Nile36N.shp');
coast_line = shaperead('d:/Egypt_EW_all/GIS/coastline_UTM_polygon.shp');
% coast_line = shaperead("D:\Egypt_EW_all\GIS\RedSea_UTM36.shp");
faults = shaperead('D:\Egypt_EW_all\GIS\Egypt faults_GeoAbdo\faults_Abdo_UTM.shp');
kharga = shaperead("D:\Egypt_EW_all\GIS\Kharga_dakhla_UTM.shp");
% farafra = shaperead("D:\Egypt_EW_all\GIS\farafra_UTM.shp");
hotSpots = readtable('D:\Egypt_EW_all\GIS\HotSpots_UTM.csv');

%
mapshow(kharga,'FaceColor','None','Linewidth',1);
% mapshow(farafra,'FaceColor','None');
% mapshow(coast_line,'FaceColor','#None','FaceAlpha',.8)
mapshow(coast_line,'FaceColor','None')

mapshow(river,'FaceColor','None' )
% mapshow(faults,'color','k','Linewidth',.2 )

plot3(hotSpots.X,hotSpots.Y, hotSpots.Z*0,'marker','pentagram','markersize',7,'color','k','linestyle','none','markerfacecolor','r');

end

%% plot rms and misfit
function plotRMS(projectPath)
misfit.file =  dir(fullfile(projectPath,['misfit.out']));
misfit.file=fullfile(misfit.file.folder,misfit.file.name);
misfit.data = readtable(misfit.file,'FileType', 'text');

rms.file =  dir(fullfile(projectPath,['rms.out']));
rms.file=fullfile(rms.file.folder,rms.file.name);
rms.data = readtable(rms.file,'FileType', 'text');

tile = tiledlayout(1,3);tile.TileSpacing = "tight";
nexttile
plot(rms.data.Var1+1,rms.data.Var2,'b');hold on
plot(rms.data.Var1+1,rms.data.Var3,'r');
title("RMS");xlabel("Iterations");ylabel("RMS");
legend({"Gravity";"MT"});

nexttile
plot(misfit.data.Var1+1,misfit.data.Var5,'b');hold on
plot(misfit.data.Var1+1,misfit.data.Var6,'r');
title("Model roughness");xlabel("Iterations");ylabel("Model roughness");
legend({"Gravity";"MT"});

nexttile
plot(misfit.data.Var1+1,misfit.data.Var7,'b');
title("MI coupling");xlabel("Iterations");ylabel("MI Coupling");
xline(misfit.data.Var1(end),'--',["Max iteration:" num2str(misfit.data.Var1(end)) ])
end

