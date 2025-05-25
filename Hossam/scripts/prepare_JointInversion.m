%% load mt and gravity files

gravity.dataFile = "D:\jif3d_JointInversion\Egypt_Grav_Mag\Data\WGM2012_20231126105429\plot_bouguer_20231126105429.txt";
gravity.dataType = "WGM2012";

%%
disp("Loading...")
gravity.dataFile = "D:\jif3d_JointInversion\All_Egypt\HH_XGM2019e.gdf";
gravity.dataType = "XGM2019e";

gravity = load_gravity_DataFile(gravity);
% shift data
fitMeshTo = "gravity"; %gravity; MT


MT.dataFile = "D:\jif3d_JointInversion\EW_profile_Egypt\prepare_inversion\data_modem_rot-15_allFreq_skip2_masked_removeSites_Zxx10_Zxy7.dat";
MT = load_MT_DataFile(MT);
%
xshift = 2*abs(min(gravity.X));
yshift = 2*abs(min(gravity.Y));
% xshift = 10*abs(min(MT.X));
% yshift = 10*abs(min(MT.Y));
% xshift = 0
% yshift = 0
gravity.xshift = xshift; gravity.yshift = yshift;
MT.xshift = xshift;MT.yshift = yshift;

for i_site = 1:size(MT.Data,2)
    MT.Data{i_site}.siteLoc(:,1:2) = [MT.X+MT.xshift,  MT.Y+MT.yshift ];
end
disp("Importing Done...✅")
%%
save('DataInfo',"MT","gravity")
disp("Saving Done...✅")

%% ploting tests latlog

scatter(gravity.long,gravity.lat,[],gravity.Data,'filled');
hold on
plot(MT.Data{1}.lon,MT.Data{1}.lat,'x','Color','k');

%% ploting tests UTM
cmap = load("D:\Software\utils\crameri colormap\CrameriColourMaps7.0.mat");
scatter(gravity.Y,gravity.X,[],gravity.Data/ 1e-5,'filled');
colormap(cmap.roma);colorbar;
hold on
plot(MT.Y,MT.X,'x','Color','k');
plot(MT.Yrot,MT.Xrot,'o','Color','k');
plotshpFiles_UTM()
% scatter(gravity.X+xshift,gravity.Y+yshift,[],gravity.Data,'filled');
%% ploting tests shifted
cmap = load("D:\Software\utils\crameri colormap\CrameriColourMaps7.0.mat");
scatter(gravity.Y+yshift,gravity.X+xshift,[],gravity.Data,'filled');
colormap(cmap.roma);colorbar;
hold on
plot(MT.Yshifted,MT.Xshifted,'x','Color','k');
%% rotate gravity data and mask data
[gravity.Yrot gravity.Xrot] = rotate_points(gravity.Y,gravity.X,MT.rotation,MT.originUTM(1),MT.originUTM(2));
gravity.Yrot = gravity.Yrot+MT.originUTM(1);
gravity.Xrot = gravity.Xrot+MT.originUTM(2);

buffer = 200000;
ind.y=find(gravity.Yrot >= min(MT.Y)-buffer & gravity.Yrot <= max(MT.Y)+buffer);
ind.x=find(gravity.Xrot >= min(MT.X)-buffer & gravity.Xrot <= max(MT.X)+buffer);
indx = intersect(ind.x,ind.y);
gravity.mask = indx;


% overwrite data coordinates with masked
[x y]= meshgrid(min(gravity.Xrot(indx)):10000:max(gravity.Xrot(indx)), ...
    min(gravity.Yrot(indx)):10000: max(gravity.Yrot(indx)));
F = scatteredInterpolant(gravity.Xrot, gravity.Yrot,gravity.Data);
d_new = F(x(:),y(:));


gravity.Y = y(:);
gravity.X = x(:);
gravity.Z = -0.3*ones(size(d_new));
gravity.Data = d_new;
gravity.dgz = ones(size(d_new)).*1e-5;


figure;nexttile
% scatter(gravity.Yrot(indx),gravity.Xrot(indx),[],gravity.Data(indx)/ 1e-5)
scatter(gravity.Y,gravity.X,[],gravity.Data/ 1e-5,'filled');hold on

scatter(y(:),x(:),[],d_new/ 1e-5,'filled');hold on
plot(MT.Y,MT.X,'x')

% nexttile
% scatter(gravity.Yrot(indx),gravity.Xrot(indx),[],gravity.Data(indx))
% scatter(y(:),x(:),[],d_new,'filled');hold on
% plot(MT.Y,MT.X,'x')
%%
plot(gravity.Y ,gravity.X,'+');
hold on
% plot(gravity.Yrot ,gravity.Xrot,'o');
plot(MT.Y,MT.X,'Marker','v','MarkerFaceColor','k','LineStyle','none')
plot(MT.Yrot,MT.Xrot,'Marker','v','MarkerFaceColor','w','LineStyle','none')

legend({"original","rotated"})
%% write grayvity data file

write_gravityDataFile(gravity);
%% write MT data file with new real coordinates
% [gravity.Y,gravity.X,~] = deg2utm( gravity.lat, gravity.long,gravity.UTMZone);
% m = [MT.sitesX 
% for i_site = 1:size(allData,2)
%     allData{i_site}.siteLoc(:,1:2) = [MT.Y+xshift MT.X+yshift];
% end

writeZ_3D('newData.dat',MT.Data,MT.header,MT.units,MT.isign)
%% prepare for gravity/MT inversion mesh

% # set mesh resolution
model_res = 15000;
xy_increaseFactor = 1.7;
z_increaseFactor = 1.3;

% #how many horizontal padding cells do we want
padx = 5;
pady = 5;
% #set the corresponding cell sizes in m
deltax = model_res;
deltay = model_res;
deltaz = 1000;

% #we do not consider topography here, to we set the top of the model to a nominal elevation of 0
startz = 0;

if fitMeshTo == "MT"
    data_x = MT.X;
    data_y = MT.Y;
else 
    data_x = gravity.X;
    data_y = gravity.Y;

end
data_spacing = sqrt((data_x(1) - data_x(2))^2 + (data_y(1) - data_y(2))^2);
fprintf('\nData spacing: %g\n', data_spacing);

% # check x and y model data extension
data_extension_x = max(data_x+xshift) - min(data_x+xshift);
data_extension_y = max(data_y+yshift) - min(data_y+yshift);


% # calculate 1 percent padding around the meassured data points
padding_x = (data_extension_x * 1)/100;
padding_y = (data_extension_y * 1)/100;

% # add padding range onto the model extention
model_extention_x = data_extension_x + (2 * padding_x); 
model_extention_y = data_extension_y + (2 * padding_y);

% #set the number of cells in each direction
nx = round(model_extention_x/model_res);if ~(iseven(nx)) ,nx = nx+1;end
ny = round(model_extention_y/model_res);if ~(iseven(ny)) ,ny = ny+1;end
nz = 18;



% #calculate the center position for the measurements
centerx = (min(data_x+xshift) + max(data_x+xshift))/2.0;
centery = (min(data_y+yshift) + max(data_y+yshift))/2.0;
% 
% % #calculate the model origin so the center of measurements is in the center of the model
% startx = centerx - nx/2 * deltax;
% starty = centery - ny/2 * deltay;
% % #calculate the layer thicknesses with depth
% #we round to full meters because
thick = round(deltaz * power(z_increaseFactor, (0:nz-1)));

% #calculate layer depth from thickness
Depth = cumsum(thick);


% #create a background density model
bg_densities =  zeros(1,nz);
bg_dens_thickness = thick;
bg_conductivities = bg_densities +.01; % 100 ohm

xpadspacing = round(deltax * (xy_increaseFactor .^ (0:(padx-1))));
ypadspacing = round(deltay * (xy_increaseFactor .^ (0:(pady-1))));
xpadwidth = sum(xpadspacing);
ypadwidth = sum(ypadspacing);

% Calculate the model origin so the center of measurements is in the center of the model
startx = centerx - nx/2 * deltax - xpadwidth;
starty = centery - ny/2 * deltay - ypadwidth;

xpadbottom = startx + cumsum(flip(xpadspacing));
ypadleft = starty + cumsum(flip(ypadspacing));

xcore = xpadbottom(end) + cumsum(ones(1, nx) * deltax);
ycore = ypadleft(end) + cumsum(ones(1, ny) * deltay);

xpadtop = xcore(end) + cumsum(xpadspacing);
ypadright = ycore(end) + cumsum(ypadspacing);
% 
Grid_X = round([xpadbottom, xcore, xpadtop]);
Grid_Y = round([ypadleft, ycore, ypadright]);
% % #create a 3D starting modell array with a priror density values
Density = zeros(nx+2*pady,ny+2*pady,nz);
padd_color.x = [zeros(size(xpadbottom)) ones(size(xcore)) zeros(size(xpadtop))];
padd_color.y = [zeros(size(ypadleft)) ones(size(ycore)) zeros(size(ypadright))];
padd_color.mesh = repmat(padd_color.y,[size(padd_color.x,2),1])
% Grid_X = [startx; Grid_X']
% Grid_Y= [starty; Grid_Y']

% % #calculate cell bounbdaries in north direction from specified number of cells and cell size
% % ## Here use grid X instead of Northing and Grid y for Easting respectivly
nx = nx+2*padx ;
ny = ny+2*pady ;
startx = centerx - nx/2 * deltax;
starty = centery - ny/2 * deltay;
Grid_X = round(startx + cumsum(ones(1, nx) * deltax));
Grid_Y = round(starty + cumsum(ones(1, ny) * deltay));
% Grid_X  = [min(Grid_X)-100000  Grid_X max(Grid_X)+100000];
% Grid_Y  = [min(Grid_Y)-100000  Grid_Y max(Grid_Y)+100000];
Density = zeros(nx, ny, nz);
% Density = zeros(nz, ny, nx);

% 

Conductivity = Density +0.01; % 100 ohm

if min(Grid_X) <0
warndlg({'detected negative cooridnaes in X(Northing) direction';'Increase the xshift value'});
elseif min(Grid_Y) <0
warndlg({'detected negative cooridnaes in Y(Easting) direction';'Increase the yshift value'});
end

disp(['nx: ' num2str(nx+2*padx) ', ny: ' num2str(ny+2*pady) ', nz: ' num2str(nz)]);
% # print thickness and depth of vertical discretization
fprintf('Thicknesses: ');
fprintf('%g ', thick);
fprintf('[m] \nDepth: ');
fprintf('%g ', Depth);
fprintf('[m]\n');
fprintf('Starting density: ');
fprintf('%g ', bg_densities);
fprintf('\n');
fprintf('Starting conuctivity: ');
fprintf('%g ', bg_conductivities);
fprintf('\n');
% plot mesh/data extends
[x y z] = meshgrid(Grid_Y,Grid_X,Depth);
% [x y z] = meshgrid([starty Grid_Y],[startx Grid_X ],[startz Depth]);
figure
% pl =scatter(gravity.Yrot+yshift,gravity.Xrot+xshift,[20],gravity.Data,'filled');view(2)
hold on
% surf(x,y,squeeze(Density(:,:,1      )));
 % surf(x,y,x);view(2);
 slice(x,y,z,Density,mean(x,"all"),[],[0])
% plot(starty-100000,startx-100000,'o','MarkerFaceColor','r','Color','r')
plot(starty,startx,'o','MarkerFaceColor','r','Color','r')
scatter(gravity.Y+yshift,gravity.X+xshift,[],gravity.Data/ 1e-5)
plot(MT.Y+MT.yshift,MT.X+MT.xshift,'x','Color','k');

set(gca,'SortMethod','childorder','ZDir','reverse')
%% write gravity mesh file

gravity.ncMeshFile = "start_grav.nc";
if isfile(gravity.ncMeshFile )
    disp('- Overwriting data file')
    delete(gravity.ncMeshFile)
end
% creat variables with dimensions
nccreate(gravity.ncMeshFile,"Northing","Dimensions", {"Northing",numel(Grid_X)}, "Format","netcdf4");
nccreate(gravity.ncMeshFile,"Northing_Origin","Format","netcdf4");

nccreate(gravity.ncMeshFile,"Easting","Dimensions", {"Easting",numel(Grid_Y)},  "Format","netcdf4");
nccreate(gravity.ncMeshFile,"Easting_Origin","Format","netcdf4");

nccreate(gravity.ncMeshFile,"Depth","Dimensions", {"Depth",numel(Depth)},   "Format","netcdf4");
nccreate(gravity.ncMeshFile,"Depth_Origin","Format","netcdf4");
% 
nccreate(gravity.ncMeshFile,"Density","Dimensions", ...
    {"Northing",numel(Grid_X),"Easting",numel(Grid_Y),"Depth",numel(Depth)}, ...
    "Format","netcdf4");


% nccreate(gravity.ncMeshFile,"Density","Dimensions", ...
    % {"Depth",numel(Depth),"Easting",numel(Grid_Y),"Northing",numel(Grid_X)}, ...
    % "Format","netcdf4");

nccreate(gravity.ncMeshFile,"bg_densities","Dimensions", {"bg_layers",numel(bg_densities)}, "Format","netcdf4");
nccreate(gravity.ncMeshFile,"bg_thicknesses","Dimensions", {"bg_layers",numel(bg_densities)}, "Format","netcdf4");

% write to variables
ncwrite(gravity.ncMeshFile,"Northing",Grid_X(:));  % x is north
ncwriteatt(gravity.ncMeshFile,"Northing","units",'m');

ncwrite(gravity.ncMeshFile,"Northing_Origin",startx);  % x is north

ncwrite(gravity.ncMeshFile,"Easting",Grid_Y(:)); % y is east
ncwriteatt(gravity.ncMeshFile,"Easting","units",'m');
ncwrite(gravity.ncMeshFile,"Easting_Origin",starty);  % x is north

ncwrite(gravity.ncMeshFile,"Depth",Depth(:));
ncwriteatt(gravity.ncMeshFile,"Depth","units",'m');
ncwrite(gravity.ncMeshFile,"Depth_Origin",startz);  % x is north

ncwrite(gravity.ncMeshFile,"Density",Density);
ncwriteatt(gravity.ncMeshFile,"Density","units",'kg/m3');


ncwrite(gravity.ncMeshFile,"bg_densities",bg_densities(:));
ncwriteatt(gravity.ncMeshFile,"bg_densities","units",'kg/m3');

ncwrite(gravity.ncMeshFile,"bg_thicknesses",bg_dens_thickness(:));
ncwriteatt(gravity.ncMeshFile,"bg_thicknesses","units",'m');

disp('Gravity mesh file Done ... ✅')

%% write MT mesh file

MT.ncMeshFile = "start_mt.nc";
if isfile(MT.ncMeshFile )
    disp('- Overwriting data file')
    delete(MT.ncMeshFile)
end
% creat variables with dimensions
nccreate(MT.ncMeshFile,"Northing","Dimensions", {"Northing",numel(Grid_X)}, "Format","netcdf4");
nccreate(MT.ncMeshFile,"Northing_Origin","Format","netcdf4");

nccreate(MT.ncMeshFile,"Easting","Dimensions", {"Easting",numel(Grid_Y)},  "Format","netcdf4");
nccreate(MT.ncMeshFile,"Easting_Origin","Format","netcdf4");

nccreate(MT.ncMeshFile,"Depth","Dimensions", {"Depth",numel(Depth)},   "Format","netcdf4");
nccreate(MT.ncMeshFile,"Depth_Origin","Format","netcdf4");

nccreate(MT.ncMeshFile,"Conductivity","Dimensions", ...
    {"Northing",numel(Grid_X),"Easting",numel(Grid_Y),"Depth",numel(Depth)}, ...
    "Format","netcdf4");

% nccreate(MT.ncMeshFile,"Conductivity","Dimensions", ...
%     {"Depth",numel(Depth),"Easting",numel(Grid_Y),"Northing",numel(Grid_X)}, ...
%     "Format","netcdf4");

nccreate(MT.ncMeshFile,"bg_conductivities","Dimensions", {"bg_layers",numel(bg_densities)}, "Format","netcdf4");
nccreate(MT.ncMeshFile,"bg_thicknesses","Dimensions", {"bg_layers",numel(bg_conductivities)}, "Format","netcdf4");

% write to variables
ncwrite(MT.ncMeshFile,"Northing",Grid_X(:));  % x is north
ncwriteatt(MT.ncMeshFile,"Northing","units",'m');

ncwrite(MT.ncMeshFile,"Northing_Origin",startx);  % x is north

ncwrite(MT.ncMeshFile,"Easting",Grid_Y(:)); % y is east
ncwriteatt(MT.ncMeshFile,"Easting","units",'m');
ncwrite(MT.ncMeshFile,"Easting_Origin",starty);  % x is north

ncwrite(MT.ncMeshFile,"Depth",Depth(:));
ncwriteatt(MT.ncMeshFile,"Depth","units",'m');
ncwrite(MT.ncMeshFile,"Depth_Origin",startz);  % x is north

ncwrite(MT.ncMeshFile,"Conductivity",Conductivity);
ncwriteatt(MT.ncMeshFile,"Conductivity","units",'S/m');


ncwrite(MT.ncMeshFile,"bg_conductivities",bg_conductivities(:));
ncwriteatt(MT.ncMeshFile,"bg_conductivities","units",'S/m');

ncwrite(MT.ncMeshFile,"bg_thicknesses",bg_dens_thickness(:));
ncwriteatt(MT.ncMeshFile,"bg_thicknesses","units",'m');

disp('MT mesh file Done ... ✅')


%% export run file
% Location of the jif3D executables on your system
execpath = "/local/jif3D_new/bin/";
% Program we want to run
progname = "grav_jointinv";
% Name of the density starting mesh
gravmesh = "start_grav.nc";
% Name of the gravity data file
gravdata = "gravdata.nc";
% Weight for standard gravity data, typically 1 in individual inversions
gravweight = 1.0;
% Weight for the magnetic data, set to zero at the moment as we do not have magnetic data yet
magweight = 0.0;
% The joint inversion can also handle other types of data but we are not interested in those right now
% so we set the corresponding weights to zero to disable them
magvectorweight = 0.0;
mtweight = 0.0;
dcweight = 0.0;
tomographyweight = 0.0;
surfacewaveweight = 0.0;
% Weight for the FTG data, set to zero at the moment as we do not have FTG data 
ftgweight = 0;
% Weight for the density regularization, higher values = smoother models
gravregularization = 1000.0;
% Maximum number of iterations
iterations = 40;

% Write all information to file
fid = fopen('run', 'w');
fprintf(fid,'%s\n',"#!/bin/bash");
fprintf(fid,'%s<<eof \n',"/local/jif3D_new/bin/jointinv --scalrelerr 0.00 --scalminerr 1e-5 --threads 32 --stochcov 1.0 --mindens -500 --maxdens 500 ");
% fprintf(fid, '%s%s<<eof \n', execpath, progname);
fprintf(fid, '%s\n', gravmesh);
fprintf(fid, '%g\n', gravweight);
fprintf(fid, '%g\n', ftgweight);
fprintf(fid, '%s\n', gravdata);
fprintf(fid, '%s\n', gravmesh);
fprintf(fid, '%g\n', magweight);
fprintf(fid, '%g\n', magvectorweight);
fprintf(fid, '%g\n', mtweight);
fprintf(fid, '%g\n', dcweight);
fprintf(fid, '%g\n', tomographyweight);
fprintf(fid, '%g\n', surfacewaveweight);
fprintf(fid, '%g\n', gravregularization);
fprintf(fid, '%d\n', iterations);
fprintf(fid, 'eof\n');
fclose(fid);

%% ploting mesh and data 
[xx yy zz] = meshgrid(Grid_Y,Grid_X,Depth);
d = permute(Density,[3,2,1]);

slice(xx-yshift,yy-xshift,zz,d,[],[2.8*10e5],[500])
set(gca,'ZDir','reverse');
plotshpFiles_UTM()
plot(gravity.Y,gravity.X,'o')
%% test xy locations (x is north)
info1 = ncdisp('D:\jif3d_JointInversion\EW_profile_Egypt\prepare_inversion\start_grav.nc');
info2 = ncdisp('start_grav.nc');

MeasPosX1 = ncread('D:\jif3d_JointInversion\EW_profile_Egypt\prepare_inversion\gravdata.nc','MeasPosX') 
MeasPosY1 = ncread('D:\jif3d_JointInversion\EW_profile_Egypt\prepare_inversion\gravdata.nc','MeasPosY') 

MeasPosY2 = ncread('gravdata.nc','MeasPosY') 
MeasPosX2 = ncread('gravdata.nc','MeasPosX') 



%%
projectPath = "R:\jointInversion\tests\m4_gravity_MT_FixXY";
iter = 66;
%%
projectPath = "R:\jointInversion\tests\m5_mt_only";
iter = 31
%%
%%
projectPath = "R:\jointInversion\All_Egypt\m1";
data_file   = 'R:\jointInversion\All_Egypt\m1\modEMdata.dat';

iter = 206
%%
gravFile.name = fullfile(projectPath,['result' num2str(iter) '.grav.inv.nc']);
mtFile.name   = fullfile(projectPath,['result' num2str(iter) '.mt.inv.nc']);

mtFile.Northing = ncread(mtFile.name,'Northing');
mtFile.Easting  = ncread(mtFile.name,'Easting');
mtFile.Depth    = ncread(mtFile.name,'Depth');
mtFile.origin   =   [ncread(mtFile.name,'Northing_Origin'), ncread(mtFile.name,'Easting_Origin'),...
    ncread(mtFile.name,'Depth_Origin')];
mtFile.resistivity = log10(1/ncread(mtFile.name,'Conductivity'));

gravFile.Northing = ncread(gravFile.name,'Northing');
gravFile.Easting  = ncread(gravFile.name,'Easting');
gravFile.Depth    = ncread(gravFile.name,'Depth');
gravFile.origin   =   [ncread(gravFile.name,'Northing_Origin'), ncread(gravFile.name,'Easting_Origin'),...
    ncread(gravFile.name,'Depth_Origin')];
gravFile.Density = ncread(gravFile.name,'Density');


[x y z] = meshgrid(mtFile.Easting,mtFile.Northing,mtFile.Depth);
disp('Done Loading')
%% Depth slice
sliceDepth = 1000;
fig = figure; fig.Position = [ 1000         563        1249         775];
tile = tiledlayout(1,2);
ax1 = nexttile;
hold on
set(ax1,'FontWeight', 'bold','FontSize',10,'Box','on','LineWidth',1.5,'ZDir','reverse');
slice(ax1, x,y,z,mtFile.resistivity,[],[],sliceDepth) % resistivity
colormap(ax1,flip(pmkmp(50,'CubicL')));
colorbar
shading flat
axis tight
ax1.Title.String = ['Resistivity model at ' num2str(sliceDepth) ' m'];

ax2 = nexttile;
hold on
set(ax2,'FontWeight', 'bold','FontSize',10,'Box','on','LineWidth',1.5,'ZDir','reverse');
slice(ax2, x,y,z,gravFile.Density,[],[],sliceDepth)
colormap(ax2,"viridis");
colorbar
ax2.Title.String = ['Density model at ' num2str(sliceDepth) ' m'];

shading flat
axis tight
Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
% ylim([2700000,3100000])
view(2)
%% load modem data

[allData,header,units,isign,originD,info] = readZ_3D(data_file,[],'Full_Impedance');
sx = []; sy = []; ch = []; sz = [];

for ind = 1 : numel(allData)
    sx = [sx; allData{ind}.siteLoc(:,1)];
    sy = [sy; allData{ind}.siteLoc(:,2)];
    sz = [sz; allData{ind}.siteLoc(:,3)];
    ch = [ch; allData{ind}.siteChar];
end
[s,i] = unique([sx,sy],'rows'); sx = sx(i); sy = sy(i); sz = sz(i); ch = ch(i,:);
origin_loc = originD';
[y0, x0,zone] = deg2utm(origin_loc(1),origin_loc(2));
sx2 = sx+x0;
sy2 = sy+y0;
%% plot stations locations
hold on
%UTM
plot3(ax1,sy2,sx2, sz-1000,'marker','v','markersize',4,'color','k','linestyle','none','markerfacecolor','k');
plot3(ax2,sy2,sx2, sz-1000,'marker','v','markersize',4,'color','k','linestyle','none','markerfacecolor','k');

%% load shapfiles
coast_line = shaperead('d:/Egypt_EW_all/GIS/coastline_UTM_polygon.shp');
mapshow(gca,coast_line,'FaceColor','#0099cc')
mapshow(ax2,coast_line,'FaceColor','#0099cc')
%%
%% slice only on stations locations
% hold on
roi = [];
% exclude station from profile
ch2 = string(strtrim(ch));
ch_ind = find(strtrim(ch2)~='soh1');

C = sortrows([sy2(ch_ind) sx2(ch_ind)],'ascend');
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

xq = .1 : .005 :numel(roi.Position(:,1));
vq = interp1(roi.Position(:,2),xq);

xq = .1 : .005 :numel(roi.Position(:,2));
vq2 = interp1(roi.Position(:,1),xq);

[xd,zd]=meshgrid(vq,squeeze(z(1,1,:)));
[yd,zd]=meshgrid(vq2,squeeze(z(1,1,:)));
%% arbityry slice plot
fig = figure; fig.Position = [ 1000         563        1249         775];
tile = tiledlayout('flow');
ax1 = nexttile;
hold on
set(ax1,'FontWeight', 'bold','FontSize',10,'Box','on','LineWidth',1.5,'ZDir','reverse','YDir','reverse');
slRes= slice(x,y,z,mtFile.resistivity,yd,xd,zd);
colormap(ax1,flip(pmkmp(50,'CubicL')));
colorbar
shading flat
% axis tight
ax1.Title.String = ['Resistivity model'];
view(0,0)

ax2 = nexttile;
hold on
set(ax2,'FontWeight', 'bold','FontSize',10,'Box','on','LineWidth',1.5,'ZDir','reverse','YDir','reverse');
slDen= slice(ax2,x,y,z,gravFile.Density,yd,xd,zd);
colormap(ax2,"viridis");
colorbar
ax2.Title.String = ['Density model' ];

shading flat
% axis tight
% Link = linkprop([ax1, ax2, ax3],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});

% ylim([2700000,3100000])
% view(3)
% sl = slice(x,y,z,gravFile.Density,yd,xd,zd);
shading flat
view(0,0)
% zlim([-10000 70000])
%% contour slice
hold on
sl = contourslice(x,y,z,mtFile.resistivity,yd,xd,zd,[0.5:0.5:4]);
%%
for i = 1:size(sl,1)
% sl(i).FaceColor = 'black';
sl(i).EdgeColor = 'white';

end
%% plot crossplot for arbitary slice
figure
tile = tiledlayout(1,2);
ax1= nexttile;
plot(mtFile.resistivity(:),gravFile.Density(:),'Marker','*','LineStyle','none')
xlabel('log10(Resistivity)')
ylabel('Density')
title('All model')


ax2 = nexttile
% plot(slRes.CData(:),slDen.CData(:),'Marker','o')
zlevels = 1 : 3 : 20;
Res = slRes.CData;
Den = slDen.CData;

Res = Res(zlevels,1:30:end); % reduce data points by 100
Den = Den(zlevels,1:30:end); % reduce data points by 100

plot(Res',Den','Marker','*','LineStyle','none')
% choose depth levels

leg = legend(num2str(squeeze(z(1,1,zlevels))/1000))
title(leg,'Depth level (km)')
title('Cross MT profile')
tile.Title.String= 'Parameter plot';
tile.Title.FontSize = 15
xlabel('log10(Resistivity)')
ylabel('Density')

%% functions



function gravity = load_gravity_DataFile(gravity)
% fixed for me
gravity.UTMZone = '36 R';

switch gravity.dataType
    case "XGM2019e"
        tmp = readtable(gravity.dataFile,'NumHeaderLines',40,'FileType','text');
    case "WGM2012"
        tmp = readtable(gravity.dataFile,'NumHeaderLines',5,'FileType','text');
end
gravity.long = tmp.Var1;
gravity.lat  = tmp.Var2;
gravity.Z    = -0.3*ones(size(tmp.Var1));    %Measurements are located 30cm above ground

gravity.dgz   =  ones(size(tmp.Var1)).*1e-5;     %set a global error of 10^-5 m/s2
gravity.Data  = tmp.Var3 * 1e-5;              %convert data from mGal to m/s2
gravity.Data  = gravity.Data-mean(gravity.Data);  %subtract mean from data to retrieve anomalies


[gravity.Y,gravity.X,~] = deg2utm( gravity.lat, gravity.long,gravity.UTMZone);

% gravity.Xshifted = gravity.X+gravity.xhsift;
% gravity.Yshifted = gravity.Y+gravity.yhsift;
end

function MT = load_MT_DataFile(MT)

[allData,header,units,isign,originD,info] = readZ_3D(MT.dataFile,[],'Full_Impedance');
sx = []; sy = []; ch = []; sz = [];
MT.header = header;MT.units =units;MT.isign =isign;
for ind = 1 : numel(allData)
    sx = [sx; allData{ind}.siteLoc(:,1)];
    sy = [sy; allData{ind}.siteLoc(:,2)];
    sz = [sz; allData{ind}.siteLoc(:,3)];
    ch = [ch; allData{ind}.siteChar];
end
[s,i] = unique([sx,sy],'rows'); sx = sx(i); sy = sy(i); sz = sz(i); ch = ch(i,:);
origin_loc = originD';
[y0, x0,zone] = deg2utm(origin_loc(1),origin_loc(2));
sx = sx+x0;
sy = sy+y0;

MT.originD = origin_loc; % latLong
MT.originUTM = [y0 x0]; %xy
MT.originUTMZone = zone; %xy
MT.sitesName = ch;
MT.sz = sz;
MT.Y = sy;
MT.X = sx;
MT.Data = allData;
MT.rotation = allData{1}.orient;
MT.T = info{1}.per;

[MT.Yrot,MT.Xrot] = rotate_points(MT.Y,MT.X,-MT.rotation,y0,x0);
MT.Yrot = MT.Yrot+y0;
MT.Xrot = MT.Xrot+x0;

% MT.Xshifted = MT.X+MT.xshift;
% MT.Yshifted = MT.Y+MT.yshift;


% for i_site = 1:size(MT.Data,2)
%     MT.Data{i_site}.siteLoc(:,1:2) = [MT.Xshifted  MT.Yshifted ];
% end

end 

function write_gravityDataFile(gravity)
% info1 = ncinfo('D:\jif3d_JointInversion\EW_profile_Egypt\prepare_inversion\gravdata.nc');
gravity.ncDtaFile = "gravdata.nc";
if isfile(gravity.ncDtaFile )
    disp('- Overwriting data file')
    delete(gravity.ncDtaFile)
end
nccreate(gravity.ncDtaFile,"MeasPosX","Dimensions",{"StationNumber",numel(gravity.X)},"Format","netcdf4");
nccreate(gravity.ncDtaFile,"MeasPosY","Dimensions",{"StationNumber",numel(gravity.X)},"Format","netcdf4");
nccreate(gravity.ncDtaFile,"MeasPosZ","Dimensions",{"StationNumber",numel(gravity.X)},"Format","netcdf4");
nccreate(gravity.ncDtaFile,"Scalar_gravity","Dimensions",{"StationNumber",numel(gravity.X)},"Format","netcdf4");
nccreate(gravity.ncDtaFile,"dGz","Dimensions",{"StationNumber",numel(gravity.X)},"Format","netcdf4");

ncwrite(gravity.ncDtaFile,"MeasPosX",gravity.X+gravity.xshift);  % x is north
ncwriteatt(gravity.ncDtaFile,"MeasPosX","units",'m')

ncwrite(gravity.ncDtaFile,"MeasPosY",gravity.Y+gravity.yshift); % y is east
ncwriteatt(gravity.ncDtaFile,"MeasPosY","units",'m')

ncwrite(gravity.ncDtaFile,"MeasPosZ",gravity.Z);
ncwriteatt(gravity.ncDtaFile,"MeasPosZ","units",'m')

ncwrite(gravity.ncDtaFile,"Scalar_gravity",gravity.Data);
ncwriteatt(gravity.ncDtaFile,"Scalar_gravity","units",'m/s2')

ncwrite(gravity.ncDtaFile,"dGz",gravity.dgz);
ncwriteatt(gravity.ncDtaFile,"dGz","units",'m/s2')

disp('Gravity data file Done ... ✅')

end

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
    mapshow(coast_line,'FaceColor','#0099cc','FaceAlpha',.8)
    mapshow(river,'FaceColor','None' )
    % mapshow(faults,'color','k','Linewidth',.2 )
    
    plot3(hotSpots.X,hotSpots.Y, hotSpots.Z*0,'marker','pentagram','markersize',7,'color','k','linestyle','none','markerfacecolor','r');

end

