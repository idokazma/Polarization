function[sigma_zz,sigma_zr] = SigmaMatrixTM;
clear all
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
lambda=1e-6;        % vacuum wavelength in meters
radius=0.05e-6;     % scatterer radius in the same physical units as lambda
rhocA=[0:20:100]*1e-6;    % array of scatterer centers in [m]
epsilonr_b=1;           % Background relative epsilon
epsilonr_s=11.4;      % Scatterer relative epsilon
mur=1;              % Relative mu. Assumed the same for background and for scatterer
%
Omegabar=[-5:0.2:5]*1e-5;   % Rotation rates normalized to optical frequency
%
Omegabars=Omegabar.*Omegabar;
Nrhoc=length(rhocA);
NOmegab=length(Omegabar);
sigma_zz=zeros(Nrhoc,NOmegab);
sigma_zr=sigma_zz;
nb=sqrt(epsilonr_b*mur);
ns=sqrt(epsilonr_s*mur);
epsilon0=8.854187817600000e-12;
mu0=1.256637061435917e-06;
mu=mu0*mur;
cLightV=1/sqrt(epsilon0*mu0);
k0=2*pi/lambda;
k0s=k0*k0;
omega=k0*cLightV;
k0nb=k0*nb;
Contrast=ns*ns-nb*nb;
Area=pi*radius*radius;
recipArea=1/Area;
Psi00I=sigI(radius,k0nb);
sigmazz00=1/(1-k0s*Contrast*0.25*1i*Psi00I*recipArea);
%sigmazz00=1.233549785*(cos(-0.1122301906)+1i*sin(-0.1122301906));  % Numerical value from Ido, for radius=0.02\mu m and n_b=1,n_s=sqrt(11.4)
%
Factor=1i*0.125*k0s*k0s;
Factor1=0.25*k0s;
p=zeros(1,Nrhoc);
p1x=p;
p1y=p;
for ir=1:Nrhoc;
    rhoc=rhocA(ir);
    [p(ir),p1x(ir),p1y(ir)]=psI(rhoc,radius,k0nb);
end
%
p=p*Factor*k0s*Contrast*sigmazz00*recipArea;
p1x=p1x*Factor1*k0s*omega*mu*Contrast*sigmazz00*recipArea;
p1y=p1y*Factor1*k0s*omega*mu*Contrast*sigmazz00*recipArea;
for ir=1:Nrhoc;
    sigma_zz(ir,:)=sigmazz00./(1-Omegabars*p(ir));
    sigma_zr(ir,:)=1i*Omegabar*p1y(ir);
end
%
% Plotting
%
rhocAN=rhocA/lambda;
nbstr=num2str(nb);
nsstr=num2str(ns);
astr=[num2str(radius/lambda),'\lambda'];
lambdastrmu=[num2str(lambda*1e6),'\mum'];
LEG=string(rhocAN);
% plotting abs(sigma_zz)
figure
for ir=1:Nrhoc;
    plot(Omegabar,abs(sigma_zz(ir,:)),'Linewidth',2)
    hold on
end
ax=gca;
ax.FontSize=14;
LGD=legend(LEG,'Fontsize',14);
title(LGD,'\rho_c [\lambda]','Fontsize',14)
xlabel('\Omega/\omega','FontSize',16)
ylabel('|\sigma_{zz}|','FontSize',16)
title(['n_b=',nbstr,', n_s=',nsstr,', r=',astr,', \lambda=',lambdastrmu],'FontSize',16)
grid on
% plotting its angle
figure
for ir=1:Nrhoc;
    plot(Omegabar,angle(sigma_zz(ir,:)),'Linewidth',2)
    hold on
end
ax=gca;
ax.FontSize=14;
LGD=legend(LEG,'Fontsize',14);
title(LGD,'\rho_c [\lambda]','Fontsize',14)
xlabel('\Omega/\omega','FontSize',16)
ylabel('\angle \sigma_{zz}','FontSize',16)
title(['n_b=',nbstr,', n_s=',nsstr,', r=',astr,', \lambda=',lambdastrmu],'FontSize',16)
grid on
%
% plotting abs(sugma_zr)
figure
for ir=1:Nrhoc;
    plot(Omegabar,abs(sigma_zr(ir,:)),'Linewidth',2)
    hold on
end
ax=gca;
ax.FontSize=14;
LGD=legend(LEG,'Fontsize',14);
title(LGD,'\rho_c [\lambda]','Fontsize',14)
xlabel('\Omega/\omega','FontSize',16)
ylabel('|\sigma_{z\rho}|','FontSize',16)
title(['n_b=',nbstr,', n_s=',nsstr,', r=',astr,', \lambda=',lambdastrmu],'FontSize',16)
grid on
% plotting its angle
figure
for ir=1:Nrhoc;
    plot(Omegabar,angle(sigma_zr(ir,:)),'Linewidth',2)
    hold on
end
ax=gca;
ax.FontSize=14;
LGD=legend(LEG,'Fontsize',14);
title(LGD,'\rho_c [\lambda]','Fontsize',14)
xlabel('\Omega/\omega','FontSize',16)
ylabel('\angle \sigma_{z\rho}','FontSize',16)
title(['n_b=',nbstr,', n_s=',nsstr,', r=',astr,', \lambda=',lambdastrmu],'FontSize',16)
grid on
end
%%
%
function [SI] = sigI(r,k0nb);
% 
% Input:
% r=  cylinder radius
% k0nb= optical wavenumber*background refraction index
% Output: 
% the integral associated with sigma00TM: Psi00TM, BUT without the i/4 !!
% 
% 
NthetaLocal=50;
NrhoLocal=50;
dt=2*pi/NthetaLocal;
dr=r/NrhoLocal;
dtdr=dt*dr;
II=zeros(NrhoLocal,NthetaLocal);
thetaLocalArray=[0:dt:2*pi-0.5*dt];
rhoLocalArray=[0.5*dr:dr:r-0.5*dr];
SelfTerm=1+i*((2/pi)*log(0.5*k0nb*dr)-1);   % self term. needs to be multiplies by rho*dr*dtheta
% Local arrays  of rho and theta, on mesh format: (local: centered at the
% cylinder center)
[RhoAL,ThetaAL]=meshgrid(rhoLocalArray,thetaLocalArray);  
XA=RhoAL.*cos(ThetaAL);    % X coordinates array. Meshgrid format
YA=RhoAL.*sin(ThetaAL);         % Y coordinates array. Meshgrid format
for ir=1:NrhoLocal;
    rL=rhoLocalArray(ir);
    for it=1:NthetaLocal;
        thetaL=thetaLocalArray(it);
        xp=rL*cos(thetaL);
        yp=rL*sin(thetaL);
        Rx=XA-xp;               % x-xprime
        Ry=YA-yp;               % y-yprime
        R=sqrt(Rx.*Rx+Ry.*Ry);
        kz=find(R==0);
        R(kz)=0.0000000000000001;
        % the integrands
        IntegrandH=besselh(0,k0nb*R);
        
%         effective_dr = sqrt(dtdr*rL/pi);
%         SelfTerm = 1+1i*((2/pi)*log(0.5*k0nb*effective_dr)-1);
        
        IntegrandH(kz)=SelfTerm;
        Itemp=rL*IntegrandH.*RhoAL;                        
        II(ir,it)=sum(Itemp(:)).*dtdr;
    end
end
SI=sum(II(:)).*dtdr;      % sigma integration (except for external factor)
end
%
function [p,p1x,p1y] = psI(rhoc,r,k0nb);
% Input:
% rhoc= distance of the cylinder center from the rotation axis
%           We assume that rhoc is in fact on the x-axis
% r=        cylinder radius [arbitrary units]
% k0nb=k_0*n_b where k_0 is the vacuum wavenumber [in 1/(units of r above)]
%           and n_b is the background refraction ibdex
% Output: 
% p,p1x,p1y= only the integrals of psi, psi1x, and psi1y
% psi=[(i*k0^4)/8]*p
% psi1x,y = [(k0^2)/4]*p1x,y
%
NthetaLocal=40;
NrhoLocal=40;
dt=2*pi/NthetaLocal;
dr=r/NrhoLocal;
dtdr=dt*dr;
II=zeros(NrhoLocal,NthetaLocal);
II1x=II;
II1y=II;
thetaLocalArray=[0:dt:2*pi-0.5*dt];
rhoLocalArray=[0.5*dr:dr:r-0.5*dr];
% Local arrays  of rho and theta, on mesh format: (local: centered at the
% cylinder center)
[RhoAL,ThetaAL]=meshgrid(rhoLocalArray,thetaLocalArray);  
XA=RhoAL.*cos(ThetaAL)+rhoc;    % X coordinates array. Meshgrid format
YA=RhoAL.*sin(ThetaAL);         % Y coordinates array. Meshgrid format
for ir=1:NrhoLocal;
    rL=rhoLocalArray(ir);
    for it=1:NthetaLocal;
        thetaL=thetaLocalArray(it);
        xp=rhoc+rL*cos(thetaL);
        yp=rL*sin(thetaL);
        Rx=XA-xp;               % x-xprime
        Ry=YA-yp;               % y-yprime
        R=sqrt(Rx.*Rx+Ry.*Ry);
        k=find(R==0);
        R(k)=0.000000000000000001;
        rhopXrho=xp*YA-yp*XA;
        rhopXrho_s=rhopXrho.*rhopXrho;
        % the integrands
        IntegrandH=besselh(0,k0nb*R);
        Itemp=rL*IntegrandH.*RhoAL;
        Integrand=Itemp.*rhopXrho_s;                        % for \psi
        II(ir,it)=sum(Integrand(:)).*dtdr;                % for \psi
        Integrand1xy=Itemp.*rhopXrho;                       % for \psi1
        Integrand1x=Integrand1xy.*Rx;       %for \psi1 - x component
        Integrand1y=Integrand1xy.*Ry;       %for \psi1 - y component
        II1x(ir,it)=sum(Integrand1x(:)).*dtdr;
        II1y(ir,it)=sum(Integrand1y(:)).*dtdr;
    end
end
p=sum(II(:)).*dtdr;      % psi (except for external factor)
p1x=sum(II1x(:)).*dtdr;   % \psi1_x except for external factor
p1y=sum(II1y(:)).*dtdr;   % \psi1_y except for external factor
end

