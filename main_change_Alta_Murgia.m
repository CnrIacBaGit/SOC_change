clear all
close all

T=12; Tinitial=2005; Tfinal=2019;
t0=Tinitial*T;
tf=Tfinal*T+T;
Nm=tf-t0;
Ny=Nm/T;

r=0.67; % 0<r<0.5 forest; 0.5<r<1 grassland; r>1 arable 
if r<=0.5
rtype=3;
else
    if r<1 
        rtype=2;
    else 
        rtype=1;
    end 
end

%1 arable, 2 grassland, 3 forest

%NPP DAL 2005 AL 2019

NPP=NPP_Alta_Murgia_2005_2019;

%NPP0=NPP(1); %2005
% gtr 12x3 arable, grassland, forest
gtrtab = zeros(12,3);
gtrtab(:,1) = [0 0 0 5/30 5/30 5/30 .5 0 0 0 0 0]';
gtrtab(:,2) = [.05 .05 .05 .05 .10 .15 .15 .10 .10 .10 .05 .05]';
gtrtab(:,3) = [.025 .025 .025 .025 .05 .05 .05 .05 .20 .20 .20 .10]';
gtr=zeros(Nm,1);
for n=1:T:Nm
    gtr(n:n+T-1)=gtrtab(:,rtype);
end


% kc 12x3 arable, (kc=0.6 per grassland, forest)


kctab = 0.6*ones(12,3);
% % Arabile Alta Murgia
 kctab(:,1) = [0.6 0.6 0.6 0.6 0.6 0.6 0.6 1 1 1 1 0.6]';
% 

kc=zeros(Nm,1);
for n=1:T:Nm
    kc(n:n+T-1)=kctab(:,rtype);
end

[Temp,rain,Tempd,pet] = weather_Alta_Murgia_2005_2019;

%[Temp,rain,Tempd,pet] = weather_Alta_Murgia_2005_2019_improved;
 
% Thornthwaite's formula for pet estimation Ld=12
I=zeros(Nm,1);
ind=Temp>0;
I(ind) = (Temp(ind)/5).^1.5;
nm = [31 28 31 30 31 30 31 31 30 31 30 31]'/30;

nmvec = ones(Nm,1); 
for n=1:T:Nm
    I(n:n+T-1) = sum(I(n:n+T-1));
    nmvec(n:n+T-1) = nm;
end


a = 6.7e-7*I.^3 - 7.7e-5*I.^2 + 1.8e-2*I + .49;
%Scommentare se si vuole usare la stima di Thornthwaite
 pet=zeros(Nm,1);
 ind=Tempd>0;
 pet(ind)= 16*nmvec(ind).*(10*Tempd(ind)./I(ind)).^a(ind);



%Alta Murgia
d=23; cly = 50;


x = 1.67*(1.85+1.60*exp(-.0786*cly));
alpha = .46/(x+1);
beta = 1/(x+1)-alpha;
delta=1-alpha-beta;
k = [10 .3 .66 .02]/T;
eta = .49;
gamma = r/(r+1);
ag=[gamma 1-gamma 0 0]';
af=[eta eta 0 1-2*eta]';


% rho
Temp0 = mean(Temp(1:T));
ka = 47.91 ./ (1+exp(106.06./(Temp+(106.06/log(46.91)-Temp0))));

[kb, Acc, M, Mb] = coeffkb(cly,d,rain,pet,Nm,Ny,T);


rho = ka.*kb.*kc;


Nb=4;
xr=30*(r-1)/r;
kc0=0.6+Nb/30*exp(xr)/(1+exp(xr));
Acc0=mean(Acc(1:T));

 if Acc0<Mb
     kb0 = 0.2 + 0.8 * (M-Acc0/(M-Mb));
 else
     kb0=1;
 end

rho0=kb0*kc0; %ka0=1

% Non standard approximation of DeltaC (dt=1)

Lam=[0 0 0 0; 0 0 0 0; alpha alpha alpha alpha; beta beta beta beta];
MM=eye(4)-Lam;
A=-MM*diag(k);
Minv=inv(MM);
At= A*Minv;
Atinv=inv(At);

DC=zeros(4,1);
DCout=DC';
tout=T;
for n=2:Ny
    for m=1:12
        tm=(n-1)*T+m;
        rhom=rho(tm);
       
        fm=rhom*A*DC+(NPP(n)/NPP(1)*gtr(tm)-rhom/T/rho0)*ag;
        if rhom<1e-16
            phi=eye(4);
        else
            phi=Atinv*( expm(rhom*At) -eye(4) )/rhom;
        end
        
        %DC=F*DC+phi*b;
        DC = DC + phi*fm; %dt=1
        %DC=DC + fm; %EE dt=1
        DCout=[DCout;DC'];
        tout=[tout;tm];
    end
end
figure(1)
hold on
Dsoc=sum(DCout,2);
meanDsoc=zeros(Ny+1);
for i=1:Ny-1
    meanD(i+1)=mean(Dsoc( (i-1)*T+1:i*T ));
end
toutyr=Tinitial+tout/12;
plot(toutyr,Dsoc,'LineWidth',1,'LineStyle','--')
hold on
plot(toutyr(1:T:end),meanD,'LineWidth',2)

title('\DeltaSoc (pet extimated by Thornthwaite’s formula)')




figure(2)
subplot(2,2,1)
plot(tout,DCout(:,1))
title('Delta dpm')
subplot(2,2,2)
plot(tout,DCout(:,2))
title('Delta rpm')
subplot(2,2,3)
plot(tout,DCout(:,3))
title('Delta bio')
subplot(2,2,4)
plot(tout,DCout(:,4))
title('Delta hum')
%         
%  
% 
% 
% %Temperature medie annuali
meanTemp=[];
meanrho=[];
for n=1:T:Nm
    meanTemp = [meanTemp; mean(Temp(n:n+T-1))];
    meanrho = [meanrho; mean(rho(n:n+T-1))];
end
% 
figure(3)
subplot(2,1,1), hold on
plot(toutyr(1:T:end),meanTemp)
title('Averaged annual tempertaures')
subplot(2,1,2), hold on
plot(toutyr(1:T:end), meanrho)
title('Averaged annual \rho(t) function')

figure(4)
hold on
plot(toutyr(1:T:end),meanTemp)
 title('Averaged annual tempertaures')




