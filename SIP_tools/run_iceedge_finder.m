clear; close all

data_dir_root=['/net/esrdata1/springer/data/'];

data_type='SSMI2';

smooth_scale = 500; % km

if data_type=='SSMI2' % SSMI NT2 25 km
   year_start=1992;
   year_stop=2007;
   dx = 25;
   max_gap = 100; 
elseif data_type=='AMSRH' % NT2 analysis 12.5 km
   year_start=2002;
   year_stop=2011;
   dx = 12.5;
   max_gap = 100; 
elseif data_type=='AMSEB' % ASI analysis 6.25 km
   grid_res=6.25;
   year_start=2002;
   year_stop=2011;
   dx = 6.25;
   max_gap = 100; 
elseif data_type=='AMSR2' % ASI analysis 3.125 km
   year_start=2012;
   year_stop=2018;
   dx=3.125;
   max_gap = 100; 
elseif data_type=='5km5d' % Mike's model output
   year_start=2009;
   year_stop=2009;
   dx=5;
   max_gap = 100; 
else 
   error('Unknown data type');
end


% smoothing parameters

T=smooth_scale/dx;
Wn=[1/T];   
[b,a]=butter(2,Wn,'low'); % lowpass

META.filtlen=T; 
META.data_type=data_type;
META.max_gap=max_gap;
META.year_start=year_start;
META.year_stop=year_stop;

META.SLAT=70;
META.SLON=0;
META.HEMI='s'; 

META.concentration=15;

SDtime=datenum(META.year_start,1,1):1:datenum(META.year_stop,12,31);
for i=1:length(SDtime)

    AE(i).SDtime=SDtime(i);
    AE(i).nsect=[]; 
    AE(i).x=[]; AE(i).y=[];
    AE(i).xsm=[]; AE(i).ysm=[];
%    AE(i).lon=[]; AE(i).lat=[];

    [IE]=find_main_ice_edge(SDtime(i),META,data_dir_root);
    if(~isempty(IE.x)); 
% raw edge
        dx   = IE.x(2:end)-IE.x(1:(end-1)); 
        dy   = IE.y(2:end)-IE.y(1:(end-1)); 
        dd   = sqrt(dx.^2+dy.^2);
        int_dd(i)=sum(dd); 

% smoothed edge
        xsm=filtfilt(b,a,IE.x);
        ysm=filtfilt(b,a,IE.y);
        dxsm=xsm(2:end)-xsm(1:(end-1));
        dysm=ysm(2:end)-ysm(1:(end-1));
        ddsm = sqrt(dxsm.^2+dysm.^2);
        int_ddsm(i)=sum(ddsm);

        AE(i).nsect=IE.nsect;
        AE(i).x=IE.x; AE(i).y=IE.y;
        AE(i).xsm=xsm; AE(i).ysm=ysm;
        %AE(i).lon=IE.lon; AE(i).lat=IE.lat;
        AE(i).rawlen=int_dd(i);
        AE(i).filtlen=int_ddsm(i);
    end
end

data_dir=[data_dir_root,'../Data/data_',META.data_type,'/'];
META.data_file = fullfile(data_dir,[META.data_type,'_iceedge','_f',num2str(META.filtlen),'_c',num2str(META.concentration),'.mat'])
save(META.data_file,'AE','META');

make_plots=1;
plot_iceedge_finder(META,make_plots);

