function plot_moa(dotsize,dotcolor,SLAT,SLON,HEMI);
%% Function to plot MOA grounding line, coastline and ice front on existing
% graphics.
%
%  If only 2 parameters passed in, it is plotted in lat/lon; if 5
%  parameters, plotted in polar stereo according to specified SLAT, SLON,
%  HEMI used by mapll.m and mapxy.m to switch between lat/lon and polar
%  stereo.
%
% plot_moa(dotsize,dotcolor,SLAT,SLON,HEMI);
%
%  Created:         29-Mar-2010 by Laurie Padman (padman@esr.org)
%  Latest update:   29-Mar-2010 (LP)
%
%=========================================================================

load ../Data/MOA_GL.mat
%load H:\Projects\MODIS_MOA\GL_files\MOA_GL.mat

if(nargin==2);  % Plotting on lat/lon
    hold on
    % Check existing longitude limits
    lonlim=get(gca,'xlim');
    dlon=0;
    if(max(lonlim)>180); % lon is (0,360);
        GL.lon(GL.lon<0)=GL.lon(GL.lon<0)+360;
        IS.lon(IS.lon<0)=IS.lon(IS.lon<0)+360;
        CO.lon(CO.lon<0)=CO.lon(CO.lon<0)+360;
    end
       
    % Add MOA boundaries
    plot(GL.lon,GL.lat,'.','color',dotcolor,'MarkerSize',dotsize); hold on
    plot(IS.lon,IS.lat,'.','color',dotcolor,'MarkerSize',dotsize); hold on
    plot(CO.lon,CO.lat,'.','color',dotcolor,'MarkerSize',dotsize); hold on
else
    [GL.x,GL.y]=mapll(GL.lat,GL.lon,SLAT,SLON,HEMI);
    [IS.x,IS.y]=mapll(IS.lat,IS.lon,SLAT,SLON,HEMI);
    [CO.x,CO.y]=mapll(CO.lat,CO.lon,SLAT,SLON,HEMI);
    
    hold on
    % Add MOA boundaries
    plot(GL.x,GL.y,'.','color',dotcolor,'MarkerSize',dotsize); hold on
    plot(IS.x,IS.y,'.','color',dotcolor,'MarkerSize',dotsize); hold on
    plot(CO.x,CO.y,'.','color',dotcolor,'MarkerSize',dotsize); hold on
end
return
