function [DATA]=read_amsre_asi(HEMI,SDtime,data_path)
% Read 6.25 km ASI (Artist) AMSR-E data

DATA.SLAT=70;DATA.SLON=0;DATA.HEMI=HEMI;
DATA.X=[]; DATA.Y=[]; DATA.seaice=[];

DV=datevec(round(SDtime));
date=[num2str(DV(1)) num2str(DV(2),'%2.2i') num2str(DV(3),'%2.2i')];
hdf_fname = ['asi-s6250-',date,'-v5i.hdf']; 

fname = fullfile(data_path,hdf_fname);

if exist(fname)
   s=hdfinfo(fullfile(data_path,hdf_fname)); 

   ice_concentration = double(hdfread(s.SDS))*100;
   xdimsize=size(ice_concentration,2);
   ydimsize=size(ice_concentration,1);

%data_description = s.SDS.Attributes.Value;

   grid_delta=6.25;
   lowleft=[-3948.71875 -3948.71875]; % these numbers come from read_amsr2.m

% Not sure this is correct.

   x=(lowleft(1))+grid_delta*(0:(xdimsize-1));
   y=(lowleft(2))+grid_delta*(0:(ydimsize-1));

   [x2d,y2d]=meshgrid(x,y);

   DATA.seaice=ice_concentration'; DATA.X=x2d';    DATA.Y=y2d';

   load('../Data/asi_625km_mask.mat','mask_asi');
   DATA.seaice(mask_asi)=NaN; 

end

