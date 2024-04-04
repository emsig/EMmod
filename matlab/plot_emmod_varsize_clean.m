% load and plots emmod output that has nx unequal ny
clear all; close all; clc; tic; 

doabs = 1; % Plot absolute value (1) or no (0)
doang = 1; % Plot phase (1) or no (0)
folder = '../emmod/'; % Folder of emmod dataset
filename = 'simplemod_11.bin'; % Filename of emmod dataset
xsize = 2000; % number of points in x-direction
ysize = 750; % number of points in y-direction
dx = 14; % sampling in x-direction
dy = 20; % smapling in y-direction
fs = 18; % Fontsize

fprintf(['Load file ',filename,'...'])
[data,xvec,yvec] = loademmod_varsize([folder,filename],xsize,ysize,dx,dy);
data = data.';
fprintf('done\n')

if doabs == 1
    figure;
    imagesc(xvec/1000,yvec/1000,log10(abs(data)));
    caxis([-16 -7]);
    colorbar
    axis image
    xlabel('offset [km]','Fontsize',fs)
    ylabel('offset [km]','Fontsize',fs)
    title(['Amplitude of ',filename],'Fontsize',fs,'interpret','none')
    set(gca,'Fontsize',fs)
end

if doang == 1
    figure;
    imagesc(xvec/1000,yvec/1000,angle(data));
    caxis([-pi pi]);
    colorbar
    axis image
    xlabel('offset [km]','Fontsize',fs)
    ylabel('offset [km]','Fontsize',fs)
    title(['Phase of ',filename],'Fontsize',fs,'interpret','none')
    set(gca,'Fontsize',fs)
end

toc