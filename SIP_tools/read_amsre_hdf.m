function [DATA]=read_amsre_hdf(HEMI,SDtime,datapath)

warning('off');
debug_flag=1;
DATA.seaice=[];
DATA.x=[];
DATA.y=[];
DATA.SLAT=70;

DV=datevec(round(SDtime));
datadir=[datapath num2str(DV(1)) '.' num2str(DV(2),'%2.2i') '.' num2str(DV(3),'%2.2i')];
D=dir([datadir '/*.hdf']);
FILE_NAME=[datadir '/' D.name];
fileinfo=dir(FILE_NAME);
e=exist(FILE_NAME);
if(isempty(fileinfo));
    if debug_flag;
      disp([FILE_NAME,' not found']);
    end
    DATA.flag=1;
    return
end
if(fileinfo.bytes==0);
    if debug_flag;
      disp([FILE_NAME,' is empty']);
    end
    DATA.flag=1;
    return
else
    DATA.flag=0;
end
file_id = hdfgd('open', FILE_NAME, 'rdonly');
if file_id == -1
   disp(['Could not open ',FILE_NAME,'...skipping'])
   return
end
DATA.infile=D.name;
DATA.fullfilename=FILE_NAME;
DATA.dir=dir;
DATA.file_id=file_id;
% Reading Data from a Data Field
if(upper(HEMI)=='N');
    GRID_NAME='NpPolarGrid12km';
    grid_id = hdfgd('attach', file_id, GRID_NAME);
    DATAFIELD_NAME='SI_12km_NH_ICECON_DAY';
    DATA.SLON=-45; DATA.HEMI=HEMI;
elseif(upper(HEMI)=='S');
    GRID_NAME='SpPolarGrid12km';
    grid_id = hdfgd('attach', file_id, GRID_NAME);
    DATAFIELD_NAME='SI_12km_SH_ICECON_DAY';
    DATA.SLON=0; DATA.HEMI=HEMI;
end
[data1, fail] = hdfgd('readfield', grid_id, DATAFIELD_NAME, [], [], []);

% Convert M-D data to 2-D data
data=double(data1);
% Transpose the data to match the map projection
data=data';

% This file contains coordinate variables that will not properly plot.
% To properly display the data, the latitude/longitude must be remapped.
[xdimsize, ydimsize, upleft, lowright, status] = hdfgd('gridinfo', grid_id);

% Detaching from the Grid Object
hdfgd('detach', grid_id);
% Closing the File
hdfgd('close', file_id);

x=(upleft(1)/1000+6.25)+12.5*(0:(xdimsize-1));
y=(upleft(2)/1000+6.25)-12.5*(0:(ydimsize-1));
[x2d,y2d]=meshgrid(x,y);

DATA.seaice=data; DATA.X=x2d;    DATA.Y=y2d;
DATA.SLAT=70;
return
