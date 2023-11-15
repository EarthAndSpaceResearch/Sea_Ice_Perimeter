function [SSMI] = read_ssmi_nt2_hig(hemi,SDtime,dx,ssmi_data_path)
% ========================================================================
% read_ssmi.m
%
% Function to read in a Southern/Northern Ocean 25-km gridded SSM/I grid
%   (NT2 algorithm), for a given Matlab time (SDtime; given in datenum 
%   format), and return a sea ice concentration (seaice) grid and polar 
%   stereographic coordinates.
%
% INPUTS:
%           hemi        char*1  hemisphere ('s' or 'n', not case-sensitive)
%           SDtime      real    Matlab time ("datenum" format)
%           dx          real    grid size (25 km for SSM/I; 12.5 km for
%                                           AMSR)
%           ssmi_data_path (char) Path to root directory for SSMI data
%                                 directories. Must include trailing '\'
% OUTPUTS in structure function SSMI:
%
%           X,Y         real    grids of polar stereo distances (km)
%           seaice      real    II x JJ pixel grid of sea ice concentration;
%                               seaice=NaN for the land mask.
%                               If hemisphere='s', (II,JJ)=(316,332).
%                               If hemisphere='n', (II,JJ)=(304,448).
%           SLAT        real    Standard latitude (70 degrees)
%           SLON        real    Stamdard longitude (0 degrees)
%
% Written by:  Laurie Padman (Earth & Space Research): padman@esr.org
%
% First created ..... 01-April-2011, modified from software provided by NSIDC (LP)
%                 from NSIDC
% USAGE:
%          [SSMI]=read_ssmi_nt2_hig(hemisphere,SDTime,dx,ssmi_data_path);
%    e.g., [SSMI]=read_ssmi_nt2_hig('s',datenum(2003,6,30),25,'j:\projects\ssmi\');
%
% General source of data and software is at:
%          http://nsidc.org/data/nsidc-0079.html
%          Also, see official data citation in ic_ssmi.m ('help ic_ssmi')
% 
% NOTES:
%        (1) Present specified array sizes are only right for 25-km grid 
%              data.
%
%==========================================================================

% Initialize the output grids for missing file
SSMI.X=[]; SSMI.Y=[]; SSMI.seaice=[];
SSMI.SLAT=[]; SSMI.SLON=[];

if(nargin<4);    % ssmi_data_path not specified
    ssmi_data_path=' ';
end

if(nargin<3);    % dx not specified
    dx=25;
    disp(' ');
    disp('dx not specified: setting dx=25 km (SSMI)');
    disp(' ');
end
hemi=lower(hemi);
hemi_uc=upper(hemi);

if(hemi=='s');
    II=316; JJ=332;
    xy_or=[-3950 -3950];
elseif(hemi=='n');
    II=304 ; JJ=448;
    xy_or=[-3850 -5350];
else
    disp('Hemisphere must be char*1 s or n. Job terminated.')
    return
end


% Parse the file name, which looks like hyearday.nt2.hig, where 
%   'h' = 's' or 'n' (hemisphere), year=year, day=day.
D=datevec(SDtime);
day=SDtime-datenum(D(1),1,1)+1;
ffname=([ssmi_data_path '/' hemi_uc num2str(D(1),'%4.4i') '/' hemi ...
       num2str(D(1),'%4.4i') num2str(day,'%3.3i') '.nt2.hig']);
ffile = dir(ffname);
%ffile=dir([ssmi_data_path '/' hemi_uc num2str(D(1),'%4.4i') '/' hemi ...
%       num2str(D(1),'%4.4i') num2str(day,'%3.3i') '.nt2.hig']);
if(isempty(ffile))
  disp(['Cannot find file: ', ffname]);
  return; 
%else
%  disp(['Found file: ', ffname]);
end
if(isempty(ffile)); return; end
fname=[ssmi_data_path '/' hemi_uc num2str(D(1),'%4.4i') '/' ffile.name];
fid = fopen(fname,'r','b');
if(fid==-1)
  disp(['Cannot find file: ', fname]);
  return; 
%else
%  disp(['Found file: ', fname]);
end

if(hemi=='s');
    L=fread(fid,[316 332],'uint16');
end
fclose(fid);
SSMI.seaice=L;
loc=find(L<0 | L>100); seaice(loc)=NaN;

% Geolocation: (find (x,y) in km
SSMI.SLAT=70; % * NOT THE SAME AS OTHER GRIDS (which are SLAT=71) *
SSMI.HDR.RE   = 6378.273;
SSMI.HDR.E2   = 0.006693883;
SSMI.HDR.E    = sqrt(SSMI.HDR.E2);
SSMI.HDR.pi   = 3.141592654;
if (hemi=='n')
    SGN = 1.0; delta = 45;
else
    SGN = -1.0; delta = 0.0;
end
SSMI.SLON=delta;

x=xy_or(1)+dx*((1:II)-0.5);
y=xy_or(2)+dx*((1:JJ)-0.5);
[X,Y]=meshgrid(x,y);
SSMI.X=X; SSMI.Y=Y;
if(SGN<0)
    SSMI.seaice=flipud((SSMI.seaice)');
else
    SSMI.seaice=(SSMI.seaice)';
end
  
return
