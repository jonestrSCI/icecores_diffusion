% Dette program (MEM_spektrum_v1) laver MEM-analyse af en dataserier (algoritme fra: Andersen 1974)


N=max(size(s));           % antal datapunkter
%p=50;                      % filterlængde



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
   
    