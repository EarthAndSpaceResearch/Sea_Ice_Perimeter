    function S=read_roms_5day(hemi,SDtime,data_path);

% Read outptut from Mike's model

   grid_file = 'so_grd.rtopo2.2.5km.try16.nc';
   x_rho = 1.e-3*ncread([data_path,grid_file],'x_rho');
   y_rho = 1.e-3*ncread([data_path,grid_file],'y_rho');
   mask_rho = ncread([data_path,grid_file],'mask_rho');
   zice = ncread([data_path,grid_file],'zice');

%Not sure if these are correct for Mike's model
   SLAT=70; SLON=0; HEMI='s';

   

   filenum = fix((SDtime-(datenum('Jan-01-2009')-2.5))/10+329);
   data_file = ['so_avg_0',num2str(filenum),'.esr.nc'];

   tt = ncread([data_path,data_file],'ocean_time');
   time = datenum('Jan-01-2009')+(tt-(9*365*86400))/86400;


   if abs(time(1)-SDtime)<=0.5 
      aice = ncread([data_path,data_file],'aice');
      S.seaice = squeeze(aice(:,:,1))*100;
      S.X=x_rho;
      S.Y=y_rho;
      S.SLAT=SLAT;
      S.SLON=SLON;
      S.HEMI=HEMI;
   %aice(~mask_rho | (zice<0))=NaN;
   elseif abs(time(2)-SDtime)<=0.5
      aice = ncread([data_path,data_file],'aice');
      S.seaice = squeeze(aice(:,:,2))*100;
      S.X=x_rho;
      S.Y=y_rho;
      S.SLAT=SLAT;
      S.SLON=SLON;
      S.HEMI=HEMI;
   else
      S.seaice= [];
   end

