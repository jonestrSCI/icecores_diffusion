clear all;
disp(' ');
disp(' ');
disp('** This Matlab script creates MEM power spectra. **');
disp(' ');
disp(' ');
bo=input('** The script will now prompt you for specifics concerning your data - press return **');
disp(' ');
disp(' ');
fid=0;
% disp('Datafilen med isotopdata skal være på formen (eksempel)');
% disp(' ');
% disp('5.65 0');
% disp('5.67 -30.24');
% disp('5.69 -32.43');
% disp(' .      .');
% disp(' .      .');
% disp(' .      .');
% disp('9.25 -29.56');
%disp(' ');
%bo=input('Tryk enter');
%disp(' ');
while fid==0;
    filnavn=input('Type in name of data file: ','s');
    disp(' ');
    %bib=input('Indtast sti til datafilen: ','s');
    bib='';
    sti=strcat(bib,filnavn);
    fid=fopen(sti,'r');
    if fid==0
        disp(' ');
        disp('Invalid file name!');
    end
end

disp(' ');
% dt=input('Indtast tykkelse af hver måling (i år): ');
% disp(' ');
% p=input('Indtast antallet af MEM-koefficienter: ');
dt=1/12;
p=30;
disp(' ');
sigma=input('Please enter estimated diffusion length (in years): ');
disp(' ');
if sigma>0
    std_noise=input('Please enter estimated standard deviation for data noise: ');
    disp(' ');
end

M=fscanf(fid,'%f',[2,100000]);  
fclose(fid);

%M=M(:,2:502);

% Dette program (MEM_spektrum) laver MEM-analyse af en dataserier (algoritme fra: Andersen 1974)

s=(M(2,2:max(size(M)))'-mean(M(2,2:max(size(M)))')); %/std(M(2,2:max(size(M)))');   % Data-serien der analyseres


N=max(size(M))-1;           % antal datapunkter
%p=50;                      % filterlængde

fn=0;                       % Tæller der styrer figur-numrene
new_try=1;
while new_try==1;       % Her begynder den løkke der gør, at man kan vælge nye værdier for p, sigma og std_noise. 

% initialisering af algoritmen

    P0=1/N*sum(s.^2);            % eq. (3)

    for t=1:N-1
        b1(t)=s(t);
        b2(t)=s(t+1);               % eq. (9)
    end

    m=1;

% Udregning af MEM-koefficienterne a(1), a(2), a(3) osv.

    nom=0;
    den=0;

    for t=1:N-m
        nom=nom+(b1(t)*b2(t));
        den=den+(b1(t)^2+b2(t)^2);
    end

    a(m)=2*nom/den;             % eq. (7)
    P(m)=P0*(1-a(m)^2);     % eq. (10)

    for m=2:p
        for t=1:m-1
            aa(t)=a(t);
        end
    
        for t=1:N-m
            b1(t)=b1(t)-aa(m-1)*b2(t);          % eq. (8)
            b2(t)=b2(t+1)-aa(m-1)*b1(t+1);
        end
   
        nom=0;
        den=0;
    
        for t=1:N-m
            nom=nom+(b1(t)*b2(t));
            den=den+(b1(t)^2+b2(t)^2);
        end
    
        a(m)=2*nom/den;                 % eq. (7)
        P(m)=P(m-1)*(1-a(m)^2);         % eq. (10)
    
        for t=1:m-1
            a(t)=aa(t)-a(m)*aa(m-t);    % eq. (5)
        end
    end   


% Nu beregnes dataseriens powerspektrum til figuren:

    M=500;                % Antal frekvenser i power-spektret
    Power=zeros(2,M);
    fmax=1/(2*dt);        % Maksimal frekvens (svingning pr aar)
    df=fmax/M;            % Frekvens-opløsningen

    f=fmax;
    for i=1:M
        Power(1,i)=f;
        real=1;
        imag=0;
        for j=1:p
            real=real-a(j)*cos(2*pi*f*dt*j);
            imag=imag+a(j)*sin(2*pi*f*dt*j);
        end
        Power(2,i)=P(p)*dt/(real^2+imag^2);
        f=f-df;
    end

% Diffusionskurve/støjkurve beregnes

    if sigma>0
        k=(0:df:fmax);

        j=max(size(k));
        DK=zeros(1,j);

        %disp(' ');
        %P_nul_svar=input('Skal P(0) defineres automatisk (tast 1) eller manuelt (tast 2): ');
        %disp(' ');
        P_nul_svar=1;
        
        if P_nul_svar==2
            P_nul=input('Indtast værdi for P(0): ');
        else
            P_nul=mean(Power(2,M-floor(M/20):M));   % Begyndelsespunkt for diffusions-kurven estimeres (ad hoc).
        end

        for i=1:j
            DK(i)=P_nul*exp(-(sigma)^2*(2*pi*k(i))^2)+(dt*std_noise^2);
        end
    end
    
    if sigma>0
    sigma_cm=sigma*10;                     % Følgende beregninger er til for at få den rigtige signaturforklaring på figuren...
    sigma_cm1=floor(sigma_cm);              % ...det er sikkert noget rigtig slam-programmering!
    sigma_cm2=round(10*(sigma_cm-sigma_cm1));
    noise0=std_noise; 
    noise1=floor(noise0);
    noise2=round(10*(noise0-noise1));
    end
    p0=p/10;
    p1=floor(p0);
    p2=round(p-p1*10);
        
    fn=fn+1;
    figure(fn);
    semilogy(Power(1,:),(Power(2,:)),'linewidth',2);
    title(filnavn);
    xlabel('Frequency [1/year]')
    ylabel('Power density [(per mille)^2*year]');
    if sigma>0
        hold on;
        semilogy(k,(DK),'r','linewidth',2);
        legend('Spectrum based on data', strcat('Modelled specturm: Sigma=0.',char(sigma_cm1+48),char(sigma_cm2+48),' year  Noise=',char(noise1+48),'.',char(noise2+48),' per mil'));
        %legend(strcat('Power m=',char(p1+48),char(p2+48)),strcat('Sigma=0.',char(sigma_cm1+48),char(sigma_cm2+48),' year  Noise=',char(noise1+48),'.',char(noise2+48),' per mil'));
        hold off;
    else
        legend(strcat('Power m=',char(p1+48),char(p2+48)));
    end
    %AXIS([0 fmax 1e-4 100])
    
    disp(' ');
    disp(' ');
    disp('Do you want to retry (with new estimates of noise and diffusion)?' );
    new_try=input('Type 1 for "yes" or 2 for "no": ');
    if new_try==1
        disp(' ');
        p=30;
        %p=input('Indtast antallet af MEM-koefficienter: ');
        %disp(' ');
        sigma=input('Please enter estimated diffusion length (in years): ');
        if sigma>0
            disp(' ');
            std_noise=input('Please enter estimated standard deviation for data noise: ');
        end
    end
end