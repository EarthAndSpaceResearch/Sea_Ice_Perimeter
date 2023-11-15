function plot_iceedge_finder(META,plots);

%warning('off');

load(META.data_file);

AE_length=length(AE);
nsect=NaN*ones(length(AE),1); int_dd=nsect; int_ddsm=nsect;
for i=1:length(AE)
   SDtime(i)=AE(i).SDtime;
   if ~isempty(AE(i).nsect);
      nsect(i) =AE(i).nsect;
      int_ddsm(i) = AE(i).filtlen;
      int_dd(i) = AE(i).rawlen;
   end
end

% limits of domain for figures
ax=[-3300 4100 -3250 4250];

fig_root_dir = '~/figures/iceedge_figures';
fig_dir = fullfile(fig_root_dir,META.data_type);

ycolors = [0  0.4470    0.7410; ...
    0.3010    0.7450    0.9330; ...
    0.4805    0.7025    0.7965; ...
    0.6600    0.6600    0.6600; ...
    0.4250    0.3860    0.4195; ...
    0.9250    0.5375    0.4240; ...
    0.9645    0.7220    0.4375; ...
    0.9290    0.6940    0.1250; ...
    0.7945    0.6770    0.3925; ...
    1.00       0.75       0.75;...
    0.8500    0.3250    0.0980; ...
    0.4680    0.4115    0.5585; ...
    0.6350    0.0780    0.1840]; ...


lcolors = [0  0.4470    0.7410; ...
    0.3010    0.7450    0.9330; ...
    0.6600    0.6600    0.6600; ...
    0.9290    0.6940    0.1250; ...
    1.00       0.75       0.75;...
    0.8500    0.3250    0.0980; ...
    0.6350    0.0780    0.1840];

figure_10_time=734403;
%figure_10_time=734767; %Sept 21, 2011

if plots
figure(9); clf; 
hold on
set(gcf,'Position',[560 323 966 854])
plot_moa(2,'k',META.SLAT,META.SLON,META.HEMI);

figure(10); clf; hold on
set(gcf,'Position',[560 323 966 854])
plot_moa(2,'k',META.SLAT,META.SLON,META.HEMI);
axis(ax);

figure(11); clf; hold on
set(gcf,'Position',[560 323 966 854])
plot_moa(2,'k',META.SLAT,META.SLON,META.HEMI);
axis(ax);

figure(12); clf; hold on
set(gcf,'Position',[560 323 966 854])
plot_moa(2,'k',META.SLAT,META.SLON,META.HEMI);
axis(ax);

figure(13); clf; hold on
set(gcf,'Position',[560 323 966 854])
plot_moa(2,'k',META.SLAT,META.SLON,META.HEMI);
axis(ax);
end


imonth=0;
imax=0;
day_max_pie=266;
pmonth=zeros(10,1);
pyear=zeros(10,1);
%for i=1:length(AE)
for i=1:AE_length
        if (AE(i).SDtime==figure_10_time & plots & 0)
           %for irad = 1:3
           %   rmax = 1000*irad; xx = [-rmax:rmax/100:rmax];
           %   plot(xx,sqrt(rmax^2-xx.^2),'g');
           %   plot(xx,-sqrt(rmax^2-xx.^2),'g');
           %end
           figure(10);
           p2=plot(AE(i).x,AE(i).y,'b','linewidth',3);
           p1=plot(AE(i).xsm,AE(i).ysm,'r-','linewidth',2);
           %l1 = legend([p1 p2],'smoothed','unsmoothed');
           set(l1,'fontsize',12,'fontname','helvetica','fontweight','demi');
           %plot(xsm,ysm,'b.','MarkerSize',2)

           figure(9);
           p2=plot(AE(i).x,AE(i).y,'b','linewidth',3);
           p1=plot(AE(i).xsm,AE(i).ysm,'r-','linewidth',2);
           l1 = legend([p1 p2],'smoothed','unsmoothed');
           set(l1,'fontsize',12,'fontname','helvetica','fontweight','demi');
           set(p1,'linewidth',8);
           set(p2,'linewidth',8);
        end
        if (rem(AE(i).SDtime-datenum('1-1-2002'),365)==day_max_pie) & META.data_type~='5km5d' & plots
           imax=imax+1;
           figure(11);
           pyear(imax)=plot(AE(i).x,AE(i).y,'color',ycolors(imax,:),'linewidth',2);

           figure(12);
           disp(['Drawing smooth contour on ',datestr(AE(i).SDtime)]);
           p3=plot(AE(i).xsm,AE(i).ysm,'b-','linewidth',2);
        end
        %if (rem(AE(i).SDtime-datenum('1-1-2002'),365)==any([104 134 165 195 226 257 287])
        if any(AE(i).SDtime-datenum('1-1-2010')==[104 134 165 195 226 257 287]) & plots
           imonth=imonth+1;
           figure(13); hold on
           %pmonth(imonth,1)=plot(AE(i).x,AE(i).y,'color',lcolors(end-imonth+1,:),'linewidth',2);
           pmonth=plot(AE(i).x,AE(i).y,'color',lcolors(end-imonth+1,:),'linewidth',2);
        end
end

if plots

if META.data_type ~='5km5d'
figure(10);
t1=title([META.data_type,' ',num2str(META.concentration),'% ice edge for ' datestr(figure_10_time)]);
x1=xlabel('X (km)');
y1=ylabel('Y (km)');
set([gca t1 x1 y1],'fontsize',14,'fontname','helvetica','fontweight','demi');
grid on
print(fullfile(fig_dir,'iceedge_image.pdf'),'-dpdf');

% same figure but zoom in
figure(9);
axis([-3250 -2250 -1000 500]);
print(fullfile(fig_dir,'iceedge_closeup_image.pdf'),'-dpdf');

figure(11);
%t1=title([META.data_type,' ',num2str(META.concentration),'% ice edge for ' datestr(SDtime(i))]);
x1=xlabel('X (km)');
y1=ylabel('Y (km)');
set([gca t1 x1 y1],'fontsize',14,'fontname','helvetica','fontweight','demi');
plot(0,0,'kx','markersize',8);
text(0+100,0,'South Pole');
if META.data_type == 'AMSRE' | META.data_type=='AMSEB'
legend(pyear,'2002','2003','2004','2005','2006','2007','2008','2009','2010','2011');
elseif META.data_type == 'AMSR2'
%legend(pyear,'2012','2013','2014','2015','2016','2017','2018');
end
print(fullfile(fig_dir,'interannual_max_iceedge.png'),'-dpng');

figure(12);
%t1=title([META.data_type,' ',num2str(META.concentration),'% ice edge for ' datestr(SDtime(i))]);
x1=xlabel('X (km)');
y1=ylabel('Y (km)');
set([gca t1 x1 y1],'fontsize',14,'fontname','helvetica','fontweight','demi');
print(fullfile(fig_dir,'smoothed_iceedge.png'),'-dpng');

xcenter=200;
ycenter=400;
RR=2.06e4/(2*pi); %radius of circle on Sept 24
xx=[-3500:100:3800];
yy = sqrt(RR^2-(xx-xcenter).^2)+ycenter;
%yy(imag(yy)~=0)=NaN;
plot(xx,yy,'g','linewidth',2);
yy = -sqrt(RR^2-(xx-xcenter).^2)+ycenter;
%yy(imag(yy)~=0)=NaN;
plot(xx,yy,'g','linewidth',2);
plot(0,0,'kx','markersize',8);
text(0+100,0,'South Pole');
plot(xcenter,ycenter,'gx','markersize',8);
text(xcenter+100,ycenter,'Center of Circle');
print(fullfile(fig_dir,'annual_iceedge.png'),'-dpng');


figure(13);
%t1=title([META.data_type,' ',num2str(META.concentration),'% ice edge for ' datestr(SDtime(i))]);
x1=xlabel('X (km)');
y1=ylabel('Y (km)');
set([gca t1 x1 y1],'fontsize',14,'fontname','helvetica','fontweight','demi');
if META.data_type == 'AMSRE' | META.data_type=='AMSEB'
l1=legend(pmonth,'Apr','May','Jun','Jul','Aug','Sep','Oct');
end
print(fullfile(fig_dir,'monthly_iceedge.png'),'-dpng');

end

figure(50); clf; orient tall; ms=2;
time=2008+(SDtime-datenum(2008,1,1))/365.25;
hs(1)=subplot(3,1,1); plot(time,nsect,'r+','MarkerSize',ms); grid on
set(gca,'ylim',[0 4]);
hs(2)=subplot(3,1,2);
plot(time,int_dd,'bo','MarkerSize',ms); hold on
plot(time,int_ddsm,'r+','MarkerSize',ms); grid on
hs(3)=subplot(3,1,3);
plot(time,int_dd./int_ddsm,'r+','MarkerSize',ms); grid on
set(gca,'ylim',[0 4]);
xlabel('Year');
for i=1:3
    set(hs(i),'xlim',[META.year_start META.year_stop+1]);
end

figure(51); clf; orient tall; ms=4;
set(gcf,'Position',[560 323 966 854])
time=2008+(SDtime-datenum(2008,1,1))/365.25;
hs(1)=subplot(2,1,1);
p1=plot(time,int_dd,'bo','MarkerSize',ms); hold on
p2=plot(time,int_ddsm,'r+','MarkerSize',ms); grid on
y1=ylabel('Length of ice edge (km)');
set(hs(1),'ylim',[0 6.e4]);
l1=legend([p1 p2],'unfiltered','filtered');
hs(2)=subplot(2,1,2);
plot(time,int_dd./int_ddsm,'k+','MarkerSize',ms); grid on
y2=ylabel('Ratio of ice edge lengths');
x1=xlabel('Year');
set(hs(2),'ylim',[1 2]);
for i=1:3
    set(hs(i),'xlim',[META.year_start META.year_stop+1]);
end

set([hs(1) hs(2)],'fontsize',12,'fontname','helvetica','fontweight','demi');
set([x1 y1 y2],'fontsize',12,'fontname','helvetica','fontweight','demi');
set(l1,'fontsize',10,'fontname','helvetica','fontweight','demi');

print(fullfile(fig_dir,'iceedge_length.pdf'),'-dpdf');

end

doy=[120:330];
tt = rem(SDtime-datenum('Jan-01-2000'),365);
good_vals = find(~isnan(int_ddsm) & ~isnan(int_dd) & tt'>=doy(1) & tt'<=doy(end));
[tt_sorted it_sorted]=sort(tt(good_vals));

good_data = int_ddsm(good_vals);
dseries_sm = good_data(it_sorted);

good_data = int_dd(good_vals);
dseries_un = good_data(it_sorted);

disp(['Total days = ',num2str(length(tt)),' Good_vals= ',num2str(length(good_vals))]);

% iterate and throw out presumed "bad" points
niter=5;
tol=0.5e4;
for i=1:niter
   [aa bb] = polyfit(tt_sorted',dseries_sm,2);
   good_vals = abs(polyval(aa,tt_sorted)'-dseries_sm)<tol*(niter+1-i);

   tt_sorted = tt_sorted(good_vals);
   dseries_sm=dseries_sm(good_vals);
   dseries_un=dseries_un(good_vals);
end
disp(['Throwing out ',num2str(length(good_data)-length(dseries_sm)),' days that are outliers']);

pie_time=tt_sorted;

figure(153);clf;
hold on;ms=2;
p_sm=plot(pie_time,dseries_sm,'ro','MarkerSize',ms); grid on

% LOESS fit
pie_smooth_loess=fLOESS(dseries_sm,0.5);
plot(pie_time,pie_smooth_loess,'k','linewidth',3);
[mm nn]=max(pie_smooth_loess);
plot(tt_sorted(nn),pie_smooth_loess(nn),'k*');
disp(['Smoothed Max perimeter = ',num2str(pie_smooth_loess(nn)),' km on day=',num2str(tt_sorted(nn))]);

%polynomial fit
[pcoef_sm bb] = polyfit(pie_time',dseries_sm,2);
%plot(doy(1:end),polyval(pcoef_sm,doy(1:end)),'y--','linewidth',4);

% unsmoothed data

p_un=plot(tt_sorted,dseries_un,'bo','MarkerSize',ms); grid on
pie_raw_loess=fLOESS(dseries_un,0.5);
disp(['Unsmoothed  perimeter = ',num2str(pie_raw_loess(nn)),' km on day=',num2str(tt_sorted(nn))]);
plot(tt_sorted,pie_raw_loess,'m','linewidth',3);
[mm nn]=max(pie_raw_loess);
plot(tt_sorted(nn),pie_raw_loess(nn),'m*');
disp(['Max perimeter = ',num2str(pie_raw_loess(nn)),' km on day=',num2str(tt_sorted(nn))]);

[pcoef_un bb] = polyfit(tt_sorted',dseries_un,5);
%plot(doy(1:end),polyval(pcoef_un,doy(1:end)),'g--','linewidth',4);

xlabel('Day of year');
ylabel('PIE length (km)');
%legend([p_sm,p_un],'smoothed','unsmoothed','location','northwest');

%title([META.data_type,' (',num2str(META.year_stop-META.year_start+1),' years)']);
text(220,5.5e4,[META.data_type,' (',num2str(META.year_start),'-',num2str(META.year_stop),')'],'fontsize',14);
text(220,5.0e4,[num2str(META.concentration),'% concentration'],'fontsize',14);
set(gca,'fontsize',12,'fontname','helvetica','fontweight','demi');
set(gca,'xlim',[100 350],'ylim',[0.e4 6.e4]);

print(fullfile(fig_dir,'PIE_annual_cycle.png'),'-dpng');

% plot as anomalies

figure(151); clf; orient tall; ms=4;
set(gcf,'Position',[560 323 966 854])
time=2008+(SDtime-datenum(2008,1,1))/365.25;

ddoy = rem(SDtime-datenum(2000,1,1),365);
good_days = ddoy>=90 & ddoy<=300;
ddoy(~good_days)=NaN;

int_dd_anomaly =  int_dd-polyval(pcoef_un,ddoy)';
int_ddsm_anomaly =  int_ddsm-polyval(pcoef_sm,ddoy)';

bad_vals = (abs(int_ddsm_anomaly)>2.e3);
int_dd_anomaly(bad_vals)=NaN;
int_ddsm_anomaly(bad_vals)=NaN;

hs(1)=subplot(2,1,1);
p2=plot(time,int_ddsm_anomaly,'r.','MarkerSize',ms); grid on
%y1=ylabel('Length of ice edge (km)');

hs(2)=subplot(2,1,2);
p1=plot(time,int_dd_anomaly,'b.','MarkerSize',ms); hold on
l1=legend([p1 p2],'unfiltered','filtered');
x1=xlabel('Year');

for i=1:2
    set(hs(i),'xlim',[META.year_start META.year_stop+1]);
end

set([hs(1) hs(2)],'fontsize',12,'fontname','helvetica','fontweight','demi');
%set([x1 y1 y2],'fontsize',12,'fontname','helvetica','fontweight','demi');
set(l1,'fontsize',10,'fontname','helvetica','fontweight','demi');

sigma_un=std(int_dd_anomaly,'omitnan');
disp(['Std. dev. of unsmoothed perimeter anomalies: ',num2str(sigma_un)])

if plots
figure(88);clf;hold on
hh=histogram(int_dd_anomaly,'Normalization','pdf');
title('Unsmoothed perimeter anomalies')
xx=0.5*(hh.BinEdges(1:end-1)+hh.BinEdges(2:end));
[cc dd] = polyfit(xx,-log(hh.Values),2);
sigma_un = 1./sqrt(cc(1));
plot(xx,1.5e-4*exp(-((xx-sqrt(cc(3)))/sigma_un).^2));

sigma_sm=std(int_ddsm_anomaly,'omitnan');
disp(['Std. dev. of smoothed perimeter anomalies: ',num2str(sigma_sm)])
figure(89);clf; hold on;
title('Smoothed perimeter anomalies')
hh=histogram(int_ddsm_anomaly,'Normalization','pdf');
xx=0.5*(hh.BinEdges(1:end-1)+hh.BinEdges(2:end));
[cc dd] = polyfit(xx,-log(hh.Values),2);
sigma_sm = 1./sqrt(cc(1));
plot(xx,4.5e-4*exp(-((xx-sqrt(cc(3)))/sigma_sm).^2));

end


