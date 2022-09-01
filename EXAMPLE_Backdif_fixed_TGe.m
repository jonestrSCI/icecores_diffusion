% Diffusion-correction code developed by Sigfus Johnsen, 
% University of Copenhagen. which uses Maximum Entropy 
% Methods (MEM) to invert an observed power density spectrum.
% Code provided in 2011 by Bo M. Vinther, University of Copenhagen
% Code adapted in 2012 by Tyler R. Jones, University of Colorado,
% to include extrema picking

% Used in the following manuscript to reconstruct seasonal
% isotopic-variability (maximum summer and minimum winter values)
% in the WAIS Divide ice core.

% Jones, T. R., Cuffey, K. M., Roberts, W. H. G., Markle, B. R., 
% Steig, E. J., Stevens, Cd18O_original.txt. M., Valdes, P. J., Fudge, T. J., 
% Sigl, M., Hughes, A. G., Morris, V., Vaughn, B. H., Garland, J., 
% Vinther, B. M., Rozmiarek, K. S., Brashear, C. A., & 
% White, J. W. C. (2022, in press). Seasonal temperatures in 
% West Antarctica during the Holocene. Nature.

% Input files include d18O_diffused.txt, to be diffusion-corrected,
% and d18O_original.txt, a hypothetical pre-diffusion signal for
% comparison purposes to the diffusion-correction and as an 
% example for understanding fitting variables (length of 
% deconvolution filter, cut-off frequency, and diffusion length)
% within a window of data.

% In ice core science, cut-off frequency is determined by 
% quantifying the frequency at which instrumental noise 
% overwhelms the climate signal, and then choosing a sufficiently 
% lower frequency to avoid this noise. The diffusion length
% is quantified by the methods in Jones et al. 2017b (see below).
% The length of deconvolution filter is usually set to 100.

% Jones, T. R., Cuffey, K. M., White, J. W. C., Steig, E. J., 
% Buizert, C., Markle, B. R., McConnell, J. R. & Sigl, M. (2017b) 
% Water isotope diffusion in the WAIS Divide ice core during 
% the Holocene and last glacial, J. Geophys. Res. Earth Surf., 
% 122.

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 14)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultTextFontSize', 14)

% Change default line widths
set(0,'defaultlinelinewidth',.5)
set(0,'defaultaxeslinewidth',.5)
set(0,'defaultpatchlinewidth',.5) 

clear all
clc
close all

%% define files
inputfil='d18O_original.txt';
outputfil='diffusion_corrected_d18O.out';

skala=0;
Dt=1/12;

%% length of deconvolution filter (a positive, whole number - i.e 100)
DataL = 100;

%% cut-off frequency (cycles/meter) or (cycles/year)
% to avoid instrumental noise or other high-frequency variations
Max_freq = 1.58;    % cycles/meter in this case

%% diffusion length (meters) or (years)
sigma = 0.237;   % meters 

%% edit carefully below
Wiener=0;

fid=fopen(inputfil,'r');
M=fscanf(fid,'%f',[2,100000]);
fclose(fid);

Int_dO18=M';  

BD_start=1;        
BD_slut=max(size(Int_dO18))-1;

Int=Int_dO18(:,1);

dt=Dt;

L=(DataL+1)*dt/2;

fDL2=floor(DataL/2);

Spektral_filter=zeros(DataL+1,fDL2);     % Spektral-opløsningsfilteret defineres for n=1 til n=floor(DataL/2)

t=-L+dt/2;
for i=1:DataL+1
    for j=1:fDL2
        Spektral_filter(i,j)=2*L*(1/(j*pi)^2)*sin(j*pi*dt/(2*L))*(sin(j*pi*(t+dt/2)/L)-sin(j*pi*(t-dt/2)/L))/((dt));
    end
    t=t+dt;
end


disp('Ice core data are extended before deconvolution...')
disp(' ');


Int_dO18_ext=zeros(max(size(Int_dO18))+DataL+2,2);       % Måleserien forlænges før filtrering

Int_dO18_ext(fDL2+2:max(size(Int_dO18))+fDL2+1,:)=Int_dO18;

MEM_extension_v2e;



disp('Deconvolution is now carried out...')
disp(' ');


recon=zeros(1+BD_slut-BD_start,2);                  % Filen som skal indeholde den diffusionskorrigerede rekonstruktion


for N=BD_start:BD_slut;                             % Rekonstruktionen udregnes

if N/100==floor(N/100)
    disp(sprintf('Number of deconvoluted samples: %d',N));
end
    
point_num=floor(N+DataL/2+2);                       % Aktuel måling der diffusionskorrigeres
    
data=zeros(DataL+3,2);

data(:,2)=Int_dO18_ext(N+1:DataL+N+3,2);
data(:,1)=Int_dO18_ext(N+1:DataL+N+3,1);


Lz=L;

sigma_t=sigma*L/Lz;                                 % Diffusionslængden omregnes fra længde til tid

step_size=dt;

Min_wavelength=1/Max_freq;

if Wiener==0
    if 2*step_size<Min_wavelength                       % Maksimalt bølgetal udregnes
        Max=floor(2*L/Min_wavelength);
    else
        Max=floor(2*L/(2*step_size));
    end
else
    if Wiener==2
        if 2*step_size<Min_wavelength                       % Maksimalt bølgetal udregnes
            Max=floor(2*L/Min_wavelength);
        else
            Max=floor(2*L/(2*step_size));
        end        
    else
        Max=floor(2*L/(2*step_size));
    end
end

a0=0;                                          

imax=DataL+2;

if Wiener>=1
    nf=2*L/Min_wavelength;
    Noise_frac=exp(-(sigma_t)^2*(pi*nf/L)^2)/Noise_freq-exp(-(sigma_t)^2*(pi*nf/L)^2);
    Noise_damp=zeros(Max,1);
    for j=1:Max                                             % Støjkorrektion beregnes for hvert bølgetal
        Noise_damp(j)=exp(-(sigma_t)^2*(pi*j/L)^2)/(Noise_frac+exp(-(sigma_t)^2*(pi*j/L)^2));  % Wiener-filtrering
    end
else
    Noise_damp=zeros(Max,1)+1;
end
    
F=zeros(imax-1,1);                                  % Tilbagediffunderingsfilter beregnes

for j=1:Max
t=-L+dt/2;    
    for i=2:imax
        if j==1
            a0=a0+(dt/(2*L))*(data(i,2));
        end
        F(i-1)=F(i-1)+Spektral_filter(i-1,j)*exp(0.5*(sigma_t)^2*(pi*j/L)^2)*Noise_damp(j); %*(dt*pi*j/(2*L))/sin(dt*pi*j/(2*L));
        t=t+dt;
    end
end
   
z=(Int_dO18_ext(point_num,1)+Int_dO18_ext(point_num-1,1))/2;    % Dybde for aktuel måling beregnes

recon(N-BD_start+1,1)=z;
recon(N-BD_start+1,2)=a0;
                                                            % Tilbagediffunderet værdi af aktuel måling beregnes
for i=2:imax                                  
    recon(N-BD_start+1,2)=recon(N-BD_start+1,2)+F(i-1)*data(i,2);
end

end

fid=fopen(outputfil,'w');
fprintf(fid,'Raw-data: %s\r\r',inputfil);
fprintf(fid,'Filter=%4.0f, ',DataL);
fprintf(fid,'MaxFreq=%3.2f, ',Max_freq);
fprintf(fid,'Scale=%1.0f (0:Tid, 1:Dybde), ',skala);
fprintf(fid,'Dt=%1.4f\r\r',Dt);
fprintf(fid,'Sigma=%1.3f ',sigma);
fprintf(fid,'Wiener=%1.0f\r\r',Wiener);
if Wiener==1
    fprintf(fid,'Max_freq-damp=%1.2f\r\r',Noise_freq);
end
fprintf(fid,' %4.3f  %4.2f  %4.2f\r',[Int_dO18(2:max(size(Int_dO18)),1),recon(:,2),Int_dO18(2:max(size(Int_dO18)),2)]');
fclose all;

Std_recon=zeros(1+BD_slut-BD_start-200,3);

for i=1:1+BD_slut-BD_start-200
    Std_recon(i,1)=mean(recon(i+99:i+100,1));
    Std_recon(i,2)=std(recon(i:i+199,2));
    Std_recon(i,3)=std(Int_dO18(i+1:i+200,2));
end
    
fid=fopen('d18O_diffused.txt','r');
isodata=fscanf(fid,'%f',[2,500]);
isodata=isodata';


%% frequency spectrum
p=30;
s=isodata(2:370,2)-mean(isodata(2:370,2));
MEM_spectrum_v1;

figure(1);

semilogy(Power(1,:),(Power(2,:)),'b');
hold on;

s=recon(:,2)-mean(recon(:,2));
MEM_spectrum_v1;

semilogy(Power(1,:),(Power(2,:)),'k');

s=Int_dO18(2:max(size(Int_dO18)),2)-mean(Int_dO18(2:max(size(Int_dO18)),2));
MEM_spectrum_v1;

semilogy(Power(1,:),(Power(2,:)),'r');
title('MEM Power Spectrum (m=30), Frequency Domain');

    xlabel('Frequency (1/m)')
    ylabel('Power density (per mil^2*m)');

legend('Original','Diffusion Corrected','Diffused');
legend('boxoff');


%% diffusion-correction
figure(2);

plot(isodata(1:end-1,1),isodata(2:370,2),'b'); hold on;
plot(recon(:,1),recon(:,2),'k');  hold on;
plot(recon(:,1),Int_dO18(2:max(size(Int_dO18)),2),'r');
title('Time Domain');
xlabel('Depth (m)');
ylabel('\delta^{18}O (per mil)');
legend('Original','Diffusion Corrected','Diffused');
legend('boxoff');


%% extrema selection
depth_tmp = recon(:,1);
iso_tmp = recon(:,2);

depth_interp = min(depth_tmp):0.005:max(depth_tmp);
iso_interp = interp1(depth_tmp,iso_tmp,depth_interp,'linear');
even_step = depth_interp(1) - depth_interp(2);


%% check interpolation
% figure(3);
% 
% plot(depth_tmp,iso_tmp,'k'); hold on
% plot(depth_interp,iso_interp,'b');


%% findpeaks algorithm
[iso_peaks,location_peaks]=findpeaks(iso_interp,'minpeakdistance',4);
depth_peaks = depth_interp(location_peaks); 

iso_flip = iso_interp(:) * (-1);
[dD_min,location_valleys]=findpeaks(iso_flip,'minpeakdistance',1);
iso_valleys = dD_min*(-1);
depth_valleys = depth_interp(location_valleys); 

figure(4);

plot(recon(:,1),recon(:,2),'k');  hold on;
plot(depth_peaks,iso_peaks,'or'), hold on;
plot(depth_valleys,iso_valleys,'ob');
title('Peak Picking');
xlabel('Depth (m)');
ylabel('\delta^{18}O (per mil)');
legend('Original','maxima','minima');
legend('boxoff');
