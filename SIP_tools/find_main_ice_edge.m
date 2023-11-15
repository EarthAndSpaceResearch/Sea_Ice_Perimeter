function [ICEEDGE]=find_main_ice_edge(SDt,META,data_dir_root);

% max_gap = 50; % size in km of biggest data gap in a single segment-
% too big and you'll jump through islands
% too small and SSMI will get kicked out;

% clear the output variables

ICEEDGE.nsect=[]; ICEEDGE.x=[]; ICEEDGE.y=[]; 

[S]=read_seaice('S',SDt,META.data_type,data_dir_root);
if(isempty(S.seaice)); 
    disp(['No file for ' datestr(SDt)]);
    return; 
end

%
cice=S.seaice; cice(cice>100)=NaN;

min_seaice_conc=META.concentration+.000001234; % make it an unusual number so we can find it in the contour array below

% plot sea ice contours
figure(1); clf; 
c=[]; h=[];
[c,h]=contour(S.X,S.Y,cice,[min_seaice_conc min_seaice_conc],'k'); 

% c([1,2],[2:end]) has the x,y coordinates of the contour line for the designated sea ice concentration
% c(1,1) is contour level = min_seaice_conc
% c(2,1) is number of vertices

if ~isempty(c)

   seg_start = find(c(1,:)==min_seaice_conc); % find indices where segments start by using the concentration contour
   seg_lengths = c(2,seg_start);

   [seg_lengths_sorted seg_start_index_sorted]=sort(seg_lengths,'descend'); % sort by length
   seg_start_sorted = seg_start(seg_start_index_sorted);


   irange=seg_start_sorted(1) + [1:seg_lengths_sorted(1)];

   ICEEDGE.x = [c(1,irange)]; 
   ICEEDGE.y = [c(2,irange)];
   ICEEDGE.nsect=1;

   POLYNYA(1).x=[];
   POLYNYA(1).y=[];

   ip=0;
   for iseg=2:length(seg_start_index_sorted)  % add additional segments
      xx1=[];yy1=[];irange=[];
      seg_length = seg_lengths_sorted(iseg);
      if seg_length>10 % min_segment % add it if it is long enough
         irange=seg_start_sorted(iseg) + [1:seg_length];
         xx1 = c(1,irange);
         yy1 = c(2,irange);
% test if it is a polynya -- closed loop
         dist0 = sqrt( (xx1(1)-xx1(end)).^2 + (yy1(1)-yy1(end)).^2);
         if dist0 < 10;
            ip=ip+1;
            POLYNYA(ip).x = xx1;
            POLYNYA(ip).y = yy1;
         else
% check distances to existing end points
            dist1 = sqrt( (xx1(1)-ICEEDGE.x(1)).^2+(yy1(1)-ICEEDGE.y(1)).^2);
            dist2 = sqrt( (xx1(end)-ICEEDGE.x(1)).^2+(yy1(end)-ICEEDGE.y(1)).^2);
            dist3 = sqrt( (xx1(1)-ICEEDGE.x(end)).^2+(yy1(1)-ICEEDGE.y(end)).^2);
            dist4 = sqrt( (xx1(end)-ICEEDGE.x(end)).^2+(yy1(end)-ICEEDGE.y(end)).^2);

            [ss tt] = min([dist1 dist2 dist3 dist4]);

            if ss < META.max_gap
               switch tt
               case 1 % new segment is to the left of existing, but needs to be reversed 
                  ICEEDGE.x = [fliplr(xx1) ICEEDGE.x];
                  ICEEDGE.y = [fliplr(yy1) ICEEDGE.y];
                  ICEEDGE.nsect = ICEEDGE.nsect+1;
               case 2 % new segment is to the left of existing
                  ICEEDGE.x = [xx1 ICEEDGE.x];
                  ICEEDGE.y = [yy1 ICEEDGE.y];
                  ICEEDGE.nsect = ICEEDGE.nsect+1;
               case 3 % new segment is to the right of existing
                  ICEEDGE.x = [ICEEDGE.x xx1];
                  ICEEDGE.y = [ICEEDGE.y yy1];
                  ICEEDGE.nsect = ICEEDGE.nsect+1;
               case 4 % new segment is to the right of existing, but needs to be reversed 
                  ICEEDGE.x = [ICEEDGE.x fliplr(xx1)];
                  ICEEDGE.y = [ICEEDGE.y fliplr(yy1)];
                  ICEEDGE.nsect = ICEEDGE.nsect+1;
               otherwise
                   error('Could not add segment for some reason')
               end % switch
           end % if ss
        end % if dist<0
   end % if c(2
end % for iseg


figure(1); hold on
plot(ICEEDGE.x,ICEEDGE.y,'r');

for ipp = 1:size(POLYNYA)
   plot(POLYNYA(ipp).x,POLYNYA(ipp).y,'g');
end
title([datestr(SDt),'; nseg=',num2str((ICEEDGE.nsect))]);
pause(0.0001)

end

return
