function [DATA]=read_amsr2_cdf(hemi,SDtime,data_dir)

DATA.seaice=[];
DATA.SLAT=70; DATA.SLON=0; DATA.HEMI=hemi;

DV=datevec(round(SDtime));

year_start = datestr(DV,10);
month_start = datestr(DV,5);
day_start = datestr(DV,7);

data_file = ['Ant_',year_start,month_start,day_start,'_res3.125_pyres.nc'];
disp(data_file);

if exist([data_dir,data_file])

  ncid = netcdf.open([data_dir,data_file]);
  ncvarid = netcdf.inqVarID(ncid,'x');

%  [xmap ymap] = mapll(lat,lon,DATA.SLAT,DATA.SLON,DATA.HEMI);
% Using 3948.71875 as reference makes the grid agree well with mapll

  x_in = 3.125.*double(netcdf.getVar(ncid,ncvarid)-1)-(3948.71875); 
  ncvarid = netcdf.inqVarID(ncid,'y');
  y_in = 3.125.*double(netcdf.getVar(ncid,ncvarid)-1)-(3948.71875);

% Not in files starting 2013/02/01
%  ncvarid = netcdf.inqVarID(ncid,'longitude');
%  lon = double(netcdf.getVar(ncid,ncvarid))/100;

%  ncvarid = netcdf.inqVarID(ncid,'latitude');
%  lat = double(netcdf.getVar(ncid,ncvarid))/100;


  ncvarid = netcdf.inqVarID(ncid,'sea_ice_concentration');
  cice_in = single(netcdf.getVar(ncid,ncvarid));
  cice_in(cice_in>=12500) = NaN;

  netcdf.close(ncid);

  cice  = cice_in/(100); % note: still goes 0->100
  [x2d y2d] = meshgrid(x_in,y_in);
  DATA.seaice=cice; DATA.X=x2d';    DATA.Y=y2d';

end

return
