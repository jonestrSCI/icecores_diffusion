% This script performs MEM-extension of data (from paper by N. Andersen (1974))
% Script to be called from main script: Back-dif_v2e.m

mean_dO18=mean(Int_dO18(2:max(size(Int_dO18)),2));
std_dO18=std(Int_dO18(2:max(size(Int_dO18)),2));

s=(Int_dO18(2:max(size(Int_dO18)),2)-mean_dO18)/std_dO18;   % Data-series to be extended

N=max(size(Int_dO18))-1;        % Size of data series af dataserien
p=DataL;                          % Size of prediction-filter

clear a;
clear P;
clear b1;
clear b2;



% initialization

P0=1/N*sum(s.^2);               % eq. (3)

for t=1:N-1
    b1(t)=s(t);
    b2(t)=s(t+1);               % eq. (9)
end

m=1;

% calculation of the coefficients (weigths)

nom=0;
den=0;

for t=1:N-m
    nom=nom+(b1(t)*b2(t));
    den=den+(b1(t)^2+b2(t)^2);
end

a(m)=2*nom/den;             % eq. (7)
P(m)=P0*(1-a(m)^2);         % eq. (10)

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
       

% Prediction-filter coefficients are now used to extend data series by p+2 samples:

s_ext=zeros(N+p+2,1);                       % Extended series
s_ext(floor(p/2)+2:N+floor(p/2)+1,1)=s;

for i=1:floor(p/2)+1
     s_ext(floor(p/2)+2-i)=sum((a').*s_ext(floor(p/2)+3-i:floor(p/2)+2+p-i));
     s_ext(N+floor(p/2)+i+1)=sum(flipud(a').*s_ext(N+floor(p/2)-p+i+1:N+floor(p/2)+i));
end

Int_dO18_ext(2:floor(p/2)+2,2)=s_ext(1:floor(p/2)+1)*std_dO18+mean_dO18;                  % Storing all data in variable used by main script
Int_dO18_ext(N+floor(p/2)+3:N+p+3,2)=s_ext(N+floor(p/2)+2:N+p+2)*std_dO18+mean_dO18;