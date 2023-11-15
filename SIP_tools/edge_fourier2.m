clear;
%close all

% smaller T value (filter) --> smaller displacements
% smaller T value --> 

data_types = {'AMSRH'}
%data_types = {'SSMI2' 'AMSRH' 'AMSEB' 'AMSR2'}

for id = 1:length(data_types);

   data_type= char(data_types(id));

   data_dir = ['/net/esrdata1/springer/data/IceEdge/data_',data_type,'/'];
   figure_root = '/net/users/springer/iceedge_figures/';
   figure_dir = [figure_root,data_type,'/'];

   switch data_type
   case 'SSMI2' % SSMI NT2 25 km
   data_name = 'SSM/I NT2';
   year_start=1992;
   year_stop=2007;
   grid_res = 25;
   T=40; min_segment=125; % min segment length -- good
   %T=80; min_segment=125; % min segment length -- good
   case 'AMSRH' % NT2 analysis 12.5 km
   data_name='AMSRE NT2';
   data_type='AMSRH'; % some older cases 
   %data_type='AMSRE'; % some older cases 
   year_start=2002;
   year_stop=2011;
   grid_res= 12.5;
   T=400; max_gap=100;
   concentration = 80;
   %T=80; min_segment=250;
   case 'AMSEB' % ASI analysis 6.25 km
   data_name='AMSRE ASI';
   year_start=2002;
   year_stop=2011;
   grid_res= 6.25;
   T=160; min_segment=500; % min segment length -- good amount of smoothing
   case 'AMSR2' % ASI analysis 3.125 km
   data_name='AMSR2 ASI';
   year_start=2012;
   year_stop=2018;
   grid_res= 3.125;
   T=320; min_segment=1000; %min segment length -- good amount of smoothing
   case '5km5d' % Mike's model output
   data_name='ROMS 5km';
   year_start=2009;
   year_stop=2009;
   grid_res= 5.;
   T=160; min_segment=625; %min segment length
   otherwise
   error('Unknown data type');
   end

   data_desc(id)={[data_name,' ',num2str(grid_res),' km']};

   filter_length=T;
   concentration=80;
%   seg_length = min_segment;

%   data_file =fullfile(data_dir,[data_type,'_iceedge_seg',num2str(seg_length),'_f',num2str(filter_length),'.mat'])
%   load(data_file);

data_dir=['/net/esrdata1/springer/data/IceEdge/data_',data_type,'/'];
data_file = fullfile(data_dir,[data_type,'_iceedge','_f',num2str(filter_length),'_c',num2str(concentration),'.mat'])
load(data_file,'AE','META');



   ntimes = length(AE);
   days = [1:ntimes];

   icd=[];

   itd=260+365*[0:round(floor((length(days)-260)/365))];
   year_count=0;
   for it = itd
      datestr(AE(it).SDtime)
      if ~isempty(AE(it).x)
         year_count=year_count+1;

         disp(['orig lengths in file: unsmoothed=',num2str(AE(it).rawlen),' smoothed=',num2str(AE(it).filtlen)]);

         x1= AE(it).x(2:end);
         x2= AE(it).x(1:end-1);
         y1= AE(it).y(2:end);
         y2= AE(it).y(1:end-1);
         ds = sqrt((x1-x2).^2+(y1-y2).^2);
         pie_length = sum(ds);
         if (any(ds<0)); warning('negative segments in unsmoothed data');end

         x1= AE(it).xsm(2:end);
         x2= AE(it).xsm(1:end-1);
         y1= AE(it).ysm(2:end);
         y2= AE(it).ysm(1:end-1);
         ds_sm = sqrt((x1-x2).^2+(y1-y2).^2);
         pie_length_sm = sum(ds_sm);

         if (any(ds_sm<0)); error('negative segments in smoothed data');end

         disp(['lengths recalculated: unsmoothed=',num2str(pie_length),' smoothed=',num2str(pie_length_sm)]);

         disp(['Unsmoothed -smoothed length, original: ',num2str(pie_length-pie_length_sm)]);

         iceedge_dist_orig = [0 cumsum(ds_sm)]; % original spacing for smoothed data

% now interpolate onto a uniformly finely spaced ds grid

         ds_equal = 1.5625/2;  % new grid spacing in m
         new_length=2^15;  % does power of 2 help with FFT?
         iceedge_dist_equal = [0:ds_equal:(new_length-1)*ds_equal];

         iceedge_xsm = interp1(iceedge_dist_orig,AE(it).xsm,iceedge_dist_equal,'linear');
         iceedge_ysm = interp1(iceedge_dist_orig,AE(it).ysm,iceedge_dist_equal,'linear');
         good_vals = ~isnan(iceedge_xsm);

% check smoothed length again
         pie_length_sm_equal = max(iceedge_dist_equal(~isnan(iceedge_xsm)));


% check regridding

%         figure(1);clf;hold on
%         plot(AE(it).xsm,AE(it).ysm);
%         plot(iceedge_xsm,iceedge_ysm,'--');


         clear displ

         for ii = 1: length(AE(it).x)
            dd = sqrt( (AE(it).xsm(ii)-AE(it).x(ii)).^2 + (AE(it).ysm(ii)-AE(it).y(ii)).^2);
% decide on the sign based on whether the unsmoothed or smoothed curve is closer to south pole
            r1 = sqrt((AE(it).x(ii)-0.2e3).^2+(AE(it).y(ii)-0.4e3).^2);
            r2 = sqrt((AE(it).xsm(ii)-0.2e3).^2+(AE(it).ysm(ii)-0.4e3).^2);
            displ(ii) = sign(r1-r2).*dd;
         end
         iceedge_displacement = interp1(iceedge_dist_orig,displ,iceedge_dist_equal,'linear');



%
%         figure(2000);clf;hold on;
%         plot(iceedge_dist_orig,displ);
%        plot(iceedge_dist_equal,iceedge_displacement,'--');

         x1= iceedge_displacement(2:end);
         x2= iceedge_displacement(1:end-1);
         ds_diff = sqrt((x1-x2).^2 + ds_equal.^2);
          
         pie_length_diff = sum(ds_diff,'omitnan');



% Hope these are the same as before the remapping; otherwise, there is a problem

         disp(['lengths recalculated: unsmoothed=',num2str(pie_length_diff),' smoothed=',num2str(pie_length_sm_equal)]);
         disp(['Unsmoothed - smoothed length, remapped: ',num2str(pie_length_diff-pie_length_sm_equal)]);

         figure(2001);hold on
         plot(pie_length-pie_length_sm, pie_length_diff-pie_length_sm_equal,'.');


% string all years together into a single "series"; for now, the years are separated by NaNs
         icd = [icd iceedge_displacement];

% store variables for Laurie

%      iceedge_dist_equal = [0:ds_equal:new_length];
%      ice_s_direct(:) = iceedge_dist_equal';
%      ice_r_direct(:,year_count) = iceedge_displacement;
%      ice_date(:,year_count) = AE(it).SDtime;
      
      end % not isempty

   end % time loop

%   d_file = fullfile([data_type,'_displacement_seg',num2str(seg_length),'_f',num2str(filter_length),'.mat']);
%   save(d_file,'ice_date','ice_s_direct','ice_r_direct');

% Parameters for welch psd estimate

% set extra values to zero before FFTs, which would be ruined by NaNs
  % icd(isnan(icd))=0;

% just take out the NaNs; no longer a power of 2
   icd_nonan=icd(~isnan(icd));

   sample_length = 1.e3; % distance over which wavenumbers are calculated
   freq_samp = sample_length/ds_equal; % sample rate per sample length; only affects labels on plots

   window_size = 2^11; % size of window used in welch averaging in index count units
   %window_size = 2^10; % size of window used in welch averaging in index count units (smaller == smoother; longer=more lowfreq res)
    
   [pw(id,:) frsamp(id,:)]=pwelch(icd_nonan,window_size,[],[],freq_samp);
   %[pw(id,:) frsamp(id,:)]=pwelch(iceedge_displacement,window_size,[],[],freq_samp);
   figure(18);hold on % uses matlab labels; for comparison with my version
   pwelch(icd_nonan,window_size,[],[],freq_samp);
   %pwelch(iceedge_displacement,window_size,[],[],freq_samp);
 
   fr_nyq = sample_length/(2*grid_res); % nyquist frequency
   fr_cutoff=max(find(frsamp(id,:)<=fr_nyq));

   pw(id,fr_cutoff:end)=NaN;

end % data_types 

figure(19);hold on
pww=plot(frsamp',log10(pw'),'linewidth',2);
xlabel(['Zonal wavenumber (cycles/(',num2str(sample_length),' km))']);
ylabel(['log_{10}(power/wavenumber)']);% (m^2/(cycles/(',num2str(sample_length),' km))']);
%ylabel(['log(Power/wavenumber) (m^2/(cycles/(',num2str(sample_length),' km))']);
%f1=plot([1 1],[-.5 3.0],'k--','linewidth',2);
%legend([pww; f1],char(data_desc{1}),char(data_desc{2}),char(data_desc{3}),char(data_desc{4}),'Filter scale');
set(gca,'xlim',[0 50]);
grid on
figure_file = fullfile(figure_root,['displacement_pspectrum_',num2str(seg_length),'_f',num2str(filter_length),'.pdf']);
print(figure_file,'-dpdf');



