% FFT over frequency to time domain to get the time domain GPR response
% 
% This script was used to create the Figure related to the third example of
% the paper "The electromagnetic response in a layered VTI medium: A new 
% look at an old problem".
clear all; clc; tic; 

% EMmod parameters
emmodfilein = '../emmod/gprloop_twointmod';
xsize = 2000; % number of points in x-direction
ysize = 2; % number of points in y-direction
dx = 0.02; % sampling in x-direction
dy = 0.02; % sampling in y-direction
freq1 = 0; % Lowest frequency minus dfreq
dfreq = 10^6; % Frequency spacing
nfreq = 2048; % Number of frequencies computed with EMmod
component = 11; % Receiver and source geometry

% Parameters concerning the GPR dataset to be generated
fc = 250*10^6; % Center Frequency for wavelet applied
newnx = 400; % Number of desired GPR receivers in x-direction
newdx = 0.02; % Desired GPR sampling in x-direction
ypos = 0; % y-offset

% Parameters for computation of theoretical arrival time of direct wave
epsr = 9; % relative electric permittivity 
zsrc = 0.0000001; % source depth
zrcv = 0.5; % receiver depth

% Parameters for computation of arrival time of reflected and refracted wave
dos = 0; % Depth Of Surface in meters
doi = 1; % Depth Of Interface in meters

% Plotting parameters
fs = 14; % Fontsize
lw = 2; % Linewidth

[data,spacevecx,spacevecy] = loademmod_varsize([emmodfilein,'_freq',num2str(freq1+dfreq),'_',num2str(component),'.bin'],xsize,ysize,dx,dy);
freqvec = linspace(0,nfreq-1,nfreq)*dfreq;
newxvec = linspace(-newnx/2,newnx/2-1,newnx)*newdx;
[xvecgrid yvecgrid] = ndgrid(spacevecx,spacevecy);
[newxvecgrid newyvecgrid] = ndgrid(newxvec,ypos);
datamat = zeros(newnx,nfreq);
fprintf('Merging data:       ')
for m = 1:nfreq
    freq1 = freq1+dfreq;
    % Loading the data
    [data,spacevecx,spacevecy] = loademmod_varsize([emmodfilein,'_freq',num2str(freq1),'_',num2str(component),'.bin'],xsize,ysize,dx,dy);
    % Add slice to cube
    datamat(:,m) = squeeze(interpn(xvecgrid,yvecgrid,data,newxvecgrid,newyvecgrid));
    fprintf('\b\b\b\b\b\b%6.2f',m/nfreq*100)
end
fprintf('\n')

fprintf('Assign wavelet...')
fwave = zeros(1,length(freqvec));
fwave(1,:) = -(freqvec/fc).^2.*exp(-(freqvec/fc).^2);
waveletmat = repmat(fwave,[newnx 1]);
datamat = datamat.*waveletmat;
fprintf('done\n')

fprintf('FFT to the time domain...')
temp = conj(flipdim(datamat,2));
temp = cat(2,zeros(newnx,1),temp);
datamat = cat(2,temp,datamat(:,2:end));
nfreq = 2*nfreq;
xxtdata = nfreq*fftshift(ifft(fftshift(datamat,2),[],2)*dfreq,2);
dt = 1/(nfreq*dfreq);
fprintf('done\n')

fprintf('Computing arrival times for direct wave...')
clight = 299792458;
vel = clight/sqrt(epsr);
arrtime = sqrt((zsrc-zrcv)^2+newxvec.^2+ypos^2)./vel;
fprintf('done\n')

fprintf('Computing arrival times for reflected wave...')
arrtimeref = sqrt((abs(zsrc-doi)+abs(zrcv-doi))^2+newxvec.^2+ypos^2)./vel;
fprintf('done\n')

fprintf('Computing arrival times for refracted wave in the air...')
% This only works if ypos = 0;
refractang = asind(vel/clight);
arrtimerefair = (abs(zsrc-dos)/cosd(refractang)+abs(zrcv-dos)/cosd(refractang))./vel;
arrtimerefair = arrtimerefair+(abs(newxvec)-abs(zsrc-dos)*tand(refractang)-abs(zrcv-dos)*tand(refractang))./clight;
fprintf('done\n')

fprintf('Prepare plotting...')
tvec = linspace(-nfreq/2,nfreq/2-1,nfreq)*dt;
gain = 1+abs((tvec*10^9).^3);
gainmat = repmat(gain,length(newxvec),1);
fprintf('done\n')

figure; 
imagesc(newxvec,tvec*10^9,(real(xxtdata).*gainmat).')
colormap('gray')
caxis([-0.5*10^13 0.5*10^13])
colorbar
hold on;
if ypos == 0
    plot(newxvec,arrtimerefair*10^9,'g','Linewidth',lw);
end
plot(newxvec,arrtime*10^9,'r','Linewidth',lw);
plot(newxvec,arrtimeref*10^9,'y','Linewidth',lw);
hold off; 
ylim([0 80]);
xlabel('Offset [m]','Fontsize',fs)
ylabel('Time [ns]','Fontsize',fs)
set(gca,'Fontsize',fs)

toc;
