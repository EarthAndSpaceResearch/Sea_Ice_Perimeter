function [S]=read_seaice(hemi,SDtime,seaice_key,data_dir_root)

if(seaice_key=='SSMI1')
    S=read_ssmi_nt1_bin(hemi,SDtime,25,ssmi_nt1_data_path);
elseif(seaice_key=='SSMI2')
    data_path=fullfile(data_dir_root,'SSMI_NT2/South/');
    S=read_ssmi_nt2_hig(hemi,SDtime,25,data_path);
elseif(seaice_key=='AMSRE')
    S=read_amsre_bin(hemi,SDtime,12.5,amsre_data_path);
elseif(seaice_key=='AMSRH'); % AMSRE in HDF format
    data_path=fullfile(data_dir_root,'AMSRE/n5eil01u.ecs.nsidc.org/AMSA/AE_SI12.003/');
    S=read_amsre_hdf(hemi,SDtime,data_path);
elseif(seaice_key=='AMSR2'); 
    data_path=fullfile(data_dir_root,'AMSR2/ftp-projects.cen.uni-hamburg.de/seaice/AMSR2/3.125km/');
    S=read_amsr2_cdf(hemi,SDtime,data_path);
elseif(seaice_key=='AMSEB'); 
    data_path=fullfile(data_dir_root,'AMSRE-ASI/ftp-projects.cen.uni-hamburg.de/seaice/AMSR-E_ASI_IceConc/hdf/s6250/');
    S=read_amsre_asi(hemi,SDtime,data_path);
elseif(seaice_key=='5km5d'); 
    data_path='/net/esrdata1/springer/CircumAnt/5km_5day/';
    S=read_roms_5day(hemi,SDtime,data_path);
else
    error('Unknown seaice data type key');
end
return
