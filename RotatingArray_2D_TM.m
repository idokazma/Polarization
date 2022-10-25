function[PVector,Efields,alphaTM]=RotatingArray_2D_TM(params, E_inc_z, Esym);
%

%
% !!!!!! e^{-i\omega t} time convention !!!!!
%
% Input:
% ======
%
% rarray - a 2D array that defines the particles LOCATIONS and parameters.
% It has Npoints rows and 2 clolumns, where Npoints is the total number of
% particles in the structure. The i-th row contains the 2 parameters of
% the i-th particle, as follows
%
% [x,y]
% 
% Here:
% x,y - the cylinder location [microns]!!!
% 
% 
%
% The cylinders polarizability (scalar) is provided by the
% function CylPolTM(aCyl,k0,epsilon_r1,epsilon_r2) below, with the obvious
% input parameters
%
%
% Output:
%
% PVector - a 1D column array, of size Npoints elements. The jth
% element gives the polarization current of the 
% j-th cylinder
%
% Efields - The same as PVector above, but it gives the fields inside each
% cylinder
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting the parameters
% lambda_micr=params.lambda;              % frequency in vacuume wavelengths [microns]
epsilon_r1=params.er_out;               % relative dielectric constant of background
% aCyl_micr=params.radius;             % Dielectric cylinders radius in microns (all cylinders are the same!)
epsilon_r2=params.er_in;            % relative Dielectric constant of cylinders
Rotan=params.OMEGA/params.omega;                 % Normalized rotation rate: Omega/omega
Rotan = -Rotan; %%% fixing sign problem. can't find where.
source_x_location_micr=0;   % source location in microns on the x axis (y=0)
x_shift_micr=0;             % shift of rotation axis in microns 
%       (structure AND source x coordinates are moved to x+x_shift)
%
lambda=params.lambda;
aCyl=params.radius;
% aCyl=params.effective_radius;
f=params.f;
omega=2*pi*f;
x_shift=x_shift_micr*1e-6;
source_location=[source_x_location_micr*1e-6+x_shift,0];
%
% rarray=array;
% rarray(:,1)=array(:,1)+x_shift_micr;
% rarray=rarray*1e-6;         % from microns to meters

rarray = [params.sca_x; params.sca_y]';

%
%
[Npoints,nc]=size(rarray);
if nc ~=2
    'coordinates/parameters are missing, only', nc, 'parameters per point'
    stop
end
%
% Start the work...
%
rhsV=zeros(Npoints,1);
PVector=zeros(Npoints);
%
%
nb=params.n_out;    % Background rotation index
k0=params.k0;
%Computing polarizability of cylinders
alphaTM=CylPolTM(params,aCyl,k0,epsilon_r1,epsilon_r2); 
InvalphaTM=1/alphaTM;

syms xsym ysym zsym k0sym
Xsym = [xsym ysym zsym];

MM=InvalphaTM*eye(Npoints,Npoints);
    for ni=1:Npoints;
        ri=rarray(ni,:);
%         Einc=IncidentE(ri,source_location,k0,Rotan,nb);
%         Einc = double(subs(Esym,[xsym,ysym,k0sym] ,[ri,k0]));
%         Einc_z = Einc(3);
        rhsV(ni)=scalar_green([params.source_loc_x;params.source_loc_y],ri',params,0);
        for nj=1:ni-1
            rj=rarray(nj,:);
            Current_Gij=GRTM(params,ri,rj,k0,Rotan,nb);
            MM(ni,nj)=-Current_Gij;
        end
        for nj=ni+1:Npoints
            rj=rarray(nj,:);
            Current_Gij=GRTM(params,ri,rj,k0,Rotan,nb);
            MM(ni,nj)=-Current_Gij;
        end
    end
MMS=sparse(MM);
rhsVS=sparse(rhsV);
PVector=MMS\rhsVS;  % The polarization current 

% Eq. 22 E = I*i/(pi*r^2*omega*delta_eps);
factor = 1i/((pi*params.radius^2)*params.omega*(epsilon_r2-epsilon_r1)*params.e0);
% Efields=PVector*InvalphaTM; % this returns Ez0
Efields = PVector*factor;
%
% Plotting
% figure;
% scatter(rarray(:,1),rarray(:,2),[],abs(PVector),'filled');
% xlabel('x [\mum]')
% ylabel('y [\mum]')
% axis equal
% colormap ('jet')
% colorbar;
% title(['|I_{pol}|, ','\Omega/\omega=',num2str(Rotan),',', ' x_s=',num2str(x_shift_micr),'\mum, ', num2str(Npoints),' Cyls'])
% figure;
% scatter(rarray(:,1),rarray(:,2),[],abs(Efields),'filled');
% xlabel('x [\mum]')
% ylabel('y [\mum]')
% axis equal
% colormap ('jet')
% colorbar;
% title(['|E|, ','\Omega/\omega=',num2str(Rotan),',', ' x_s=',num2str(x_shift_micr),'\mum, ', num2str(Npoints),' Cyls'])
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[IncFields]=IncidentE(r0,rpi,k0,Rotan,nindex) % 
%   e^{-i\omega t} time convention
%
%rpi= Source location
factor=1/(1i*k0*3e8*1.256637061435917e-06);
IncFields=factor*GRTM(r0,rpi,k0,Rotan,nindex);    % Incident E at ro due to a 
%   line current of one unit located at rpi

end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[G]=GRTM(params, ro,rp,k0,Rotan,nindex)
% a function that creates the approximate Green's function of a 2D rotating medium in the TM pol.      
%
% Input:
% ro= 2-elements ROW vector, defines the x,y components of observation point
% rp= the same, for source location
% k0=vaccum wavenumber
% Rotan = normalized rotation rate: Omega/omega
% nindex = background refrection index
%
% Output:
% G = the response to a z-directed unit electric current (1 Amper) in an
% infinitesimal thin line. G is based on the approximate Green's function
% computation in rotating medium (Eq (3.19) in my Green's function paper,
% multiplied by i\omega\mu0)
%
mu0=1.256637061435917e-06;    % mu0
mu0 = params.mu0;
%c=3e8;
k0 = params.k0;
r=norm(ro-rp);
krn=k0*r*nindex;
k0s=k0*k0;
omega=k0*params.c;
Gst=-0.25*omega*mu0*besselh(0,1,krn);   % Stationary medium G
z_d_rp_cross_ro=diff(ro.*(flip(rp)));
RexpFactor=k0s*Rotan*z_d_rp_cross_ro;
gr=exp(1i*RexpFactor);
G=Gst*gr;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[alphaTM]=CylPolTM(params, aCyl,k0,epsilon_r1,epsilon_r2);
%
% A funcrion that computes the electric polarizability of thin dielectic
% cylinder in TM illomination (E field along the cylinder axis)
omega_mu0=params.k0*params.c*params.mu0; %k0*3e8*1.256637061435917e-06;

n1=params.n_out;
n2=params.n_in;
k0n1a=k0*n1*aCyl;
k0n2a=k0*n2*aCyl;
j0k0n2a=besselj(0,k0n2a);
j1k0n2a=besselj(1,k0n2a);
Dnum=-besselj(1,k0n1a)*j0k0n2a*n1+j1k0n2a*besselj(0,k0n1a)*n2;
Den=-besselh(1,1,k0n1a)*j0k0n2a*n1+j1k0n2a*besselh(0,1,k0n1a)*n2;
minus_b0TM=Dnum/Den;
alphaTM=4*minus_b0TM/omega_mu0;  
end
    
    
    
    
    
    