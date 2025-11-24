function[PVector,Efields,alphaTE]=RotatingArray_2D_TE_v2(params,field_current_mat);
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
epsilon_r1=params.er_out;               % relative dielectric constant of background
epsilon_r2=params.er_in;            % relative Dielectric constant of cylinders
rarray = [params.sca_x; params.sca_y]'; %scatterers array;
source_location=[params.source_loc_x,params.source_loc_y];
Rotan = params.OMEGA/params.omega; %Rotation ratio (OMEGA/omega)


[Npoints,nc]=size(rarray);
if nc ~=2
    'coordinates/parameters are missing, only', nc, 'parameters per point'
    stop
end
%
% Start the work...
%
PVector=zeros(Npoints);

%Computing polarizability of cylinders
alphaTE=CylPolTE(params); 
InvalphaTE=1/alphaTE;

GGxx = zeros(length(params.sca_x));
GGyy = zeros(length(params.sca_x));
GGyx = zeros(length(params.sca_x));
GGxy = zeros(length(params.sca_x));

for ii=1:length(params.sca_x) 
    % go over all scatterers and calculate for each two things: 
    % 1 - the incident field from the source.
    % 2 - the impact of the other scatterers radiation (green's function).
    % (the solution we are looking for is the amplitude and phase to
    % multiply this greens function.
    % we do that to solve Eq. 19 (and Eq. 18);
    [incE_x(ii), incE_y(ii)]=IncidentE(params, [params.sca_x(ii),params.sca_y(ii)],source_location);
%      incE_x(ii) = field_current_mat{3};
%      incE_y(ii) = field_current_mat{4};
     [GGxx(:,ii),GGxy(:,ii),GGyx(:,ii),GGyy(:,ii)] = dyiadic_green_matrix([params.sca_x(ii);params.sca_y(ii)],[params.sca_x;params.sca_y],params.n_out,params.k0,params.c,params.OMEGA,params); % only Hankel and Exponent part

%     GG(:,ii) = scalar_green([params.sca_x(ii);params.sca_y(ii)],[params.sca_x;params.sca_y],params.n_out,params.k0,params.c,params.OMEGA); % only Hankel and Exponent part
%     GG(:,ii) = -1 * 1i * params.omega *params.mu0 * GG(:,ii); 
%     GG(ii,ii) = InvalphaTM; % the "self term" -> I_pol / alpha_TM = E_local;
      GGxx(ii,ii) = InvalphaTE;
      GGyy(ii,ii) = InvalphaTE;
      GGyx(ii,ii) = 0;
      GGxy(ii,ii) = 0;

end
     incE_x = field_current_mat{3}.';
     incE_y = field_current_mat{4}.';

GG = [GGxx,GGxy;GGyx,GGyy];
incE = [incE_x,incE_y];
PVector = linsolve(GG,incE.');
Px = PVector(1:length((params.sca_x)));
Py = PVector(length((params.sca_x))+1:end);

% Eq. 22 E = I*i/(pi*r^2*omega*delta_eps);
factor = 1i/((pi*params.radius^2)*params.omega*(epsilon_r2-epsilon_r1)*params.e0);
% Efields=PVector*InvalphaTM; % this returns Ez0
Efields = PVector*factor;
Ex = Efields(1:length((params.sca_x)));
Ey = Efields(length((params.sca_x))+1:end);


%
% Plotting
if (params.toplot == 1)
    figure;
    scatter(rarray(:,1),rarray(:,2),[],sqrt(abs(Px).^2 + abs(Py).^2),'filled');
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    axis equal
    colormap ('jet')
    colorbar;
    title(['|I_{pol}|, ','\Omega/\omega=',num2str(Rotan),',', num2str(Npoints),' Cyls'])
    figure;
    scatter(rarray(:,1),rarray(:,2),[],sqrt(abs(Ex).^2 + abs(Ey).^2),'filled');
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    axis equal
    colormap ('jet')
    colorbar;
    title(['|E|, ','\Omega/\omega=',num2str(Rotan),',', num2str(Npoints),' Cyls'])
        caxis ([0,3e8]);

    figure;
    subplot(1,3,1);
    quiver(rarray(:,1),rarray(:,2),real(incE_x).',real(incE_y).'); title ('Hit real');
    hold on;
    quiver(rarray(:,1),rarray(:,2),imag(incE_x).',imag(incE_y).'); title ('Hit imag');
    subplot(1,3,3);
    quiver(rarray(:,1),rarray(:,2),abs(incE_x).',abs(incE_y).'); title ('Hit abs');
    
    figure;
    subplot(1,3,1);
    quiver(rarray(:,1),rarray(:,2),real(Px),real(Py)); title ('real');
    subplot(1,3,2);
    quiver(rarray(:,1),rarray(:,2),imag(Px),imag(Py)); title ('imag');
    subplot(1,3,3);
    quiver(rarray(:,1),rarray(:,2),abs(Px),abs(Py)); title ('abs');
    
    
     figure;
    scatter(rarray(:,1),rarray(:,2),[],sqrt(imag(incE_x).^2 + imag(incE_y).^2),'filled');
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    axis equal
    colormap ('jet')
    colorbar;
    title(['imag power hit ','\Omega/\omega=',num2str(Rotan),',', num2str(Npoints),' Cyls'])
    caxis ([0,3e8]);
    
         figure;
    scatter(rarray(:,1),rarray(:,2),[],sqrt(real(incE_x).^2 + real(incE_y).^2),'filled');
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    axis equal
    colormap ('jet')
    colorbar;
    title(['real power hit ','\Omega/\omega=',num2str(Rotan),',', num2str(Npoints),' Cyls'])
    caxis ([0,3e8]);
    
      
     figure;
    scatter(rarray(:,1),rarray(:,2),[],sqrt(imag(Ex).^2 + imag(Ey).^2),'filled');
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    axis equal
    colormap ('jet')
    colorbar;
    title(['|E|, ','\Omega/\omega=',num2str(Rotan),',', num2str(Npoints),' Cyls'])
%      caxis ([0,max(sqrt(imag(Ex).^2 + imag(Ey).^2))/8]);
    
         figure;
    scatter(rarray(:,1),rarray(:,2),[],sqrt(real(Ex).^2 + real(Ey).^2),'filled');
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    axis equal
    colormap ('jet')
    colorbar;
    title(['|E|, ','\Omega/\omega=',num2str(Rotan),',', num2str(Npoints),' Cyls'])
%      caxis ([0,max(sqrt(real(Ex).^2 + real(Ey).^2))/2]);



     figure;
    scatter(rarray(:,1),rarray(:,2),[],abs(incE_x),'filled');
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    axis equal
    colormap ('jet')
    colorbar;
    title(['abs x hit ','\Omega/\omega=',num2str(Rotan),',', num2str(Npoints),' Cyls'])
    
         figure;
    scatter(rarray(:,1),rarray(:,2),[],abs(incE_y),'filled');
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    axis equal
    colormap ('jet')
    colorbar;
    title(['abs y hit ','\Omega/\omega=',num2str(Rotan),',', num2str(Npoints),' Cyls'])
    
    
         figure;
    scatter(rarray(:,1),rarray(:,2),[],abs(Ex),'filled');
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    axis equal
    colormap ('jet')
    colorbar;
    title(['abs x SOL ','\Omega/\omega=',num2str(Rotan),',', num2str(Npoints),' Cyls'])
    
         figure;
    scatter(rarray(:,1),rarray(:,2),[],abs(Ey),'filled');
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    axis equal
    colormap ('jet')
    colorbar;
    title(['abs y SOL ','\Omega/\omega=',num2str(Rotan),',', num2str(Npoints),' Cyls'])
   
    
    
end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[Inc_Ex, Inc_Ey]=IncidentE(params,r0,rpi) % 
%   e^{-i\omega t} time convention
%
%rpi= Source location
factor=(1i*params.k0*params.c*params.mu0)*1;  % 1 ampere;
direction = pi/2; %0 = x axis. [0,2*pi];
phase_diff = 0;
I_x = 1*cos(direction);
I_y = 1*sin(direction);

I = [I_x;I_y];

% IncFields=factor*scalar_green(rpi',r0',params.n_out,params.k0,params.c,params.OMEGA);    % Incident E at ro due to a 
%   line current of one unit located at rpi




[Gxx,Gxy,Gyx,Gyy]=dyiadic_green_matrix(rpi',r0',params.n_out,params.k0,params.c,params.OMEGA, params);    % Incident E at ro due to a 
Inc_Ex = [Gxx,Gxy]*I;
Inc_Ey = [Gyx,Gyy]*I;

end


function[alphaTE]=CylPolTE(params);
%
% % % A funcrion that computes the electric polarizability of thin dielectic
% % % cylinder in TM illomination (E field along the cylinder axis)
% % omega_mu0=params.k0*params.c*params.mu0; %k0*3e8*1.256637061435917e-06;
% % 
% % n1=params.n_out;
% % n2=params.n_in;
% % k0n1a=params.k0*n1*params.radius;
% % k0n2a=params.k0*n2*params.radius;
% % j0k0n2a=besselj(0,k0n2a);
% % j1k0n2a=besselj(1,k0n2a);
% % Dnum=-besselj(1,k0n1a)*j0k0n2a*n1+j1k0n2a*besselj(0,k0n1a)*n2;
% % Den=-besselh(1,1,k0n1a)*j0k0n2a*n1+j1k0n2a*besselh(0,1,k0n1a)*n2;
% % minus_b0TM=Dnum/Den;
% % alphaTM=4*minus_b0TM/omega_mu0;  




n_1 = params.n_out;
n_2 = params.n_in;

k0 = params.k0;
radius = params.radius;

etha_1 = sqrt(params.mr_out/params.er_out);
etha_2 = sqrt(params.mr_in/params.er_in);

m=0:3;
tic
%% b_te calc
b_m_te_up_A = etha_1*(1/2*(besselj(m-1,k0*n_1*radius)-besselj(m+1,k0*n_1*radius))).*besselj(m,k0*n_2*radius);
b_m_te_up_B = etha_2*(1/2*(besselj(m-1,k0*n_2*radius)-besselj(m+1,k0*n_2*radius))).*besselj(m,k0*n_1*radius);

b_m_te_down_C = etha_1*(1/2*(besselh(m-1,1,k0*n_1*radius)-besselh(m+1,1,k0*n_1*radius))).*besselj(m,k0*n_2*radius);
b_m_te_down_D = etha_2*(1/2*(besselj(m-1,k0*n_2*radius)-besselj(m+1,k0*n_2*radius))).*besselh(m,1,k0*n_1*radius);

b_TE = ((-1)*(-1i).^m).*(b_m_te_up_A-b_m_te_up_B)./(b_m_te_down_C-b_m_te_down_D);

%% c_te calc

% c_TE = b_TE.*besselh(m,1,k0*n_1*radius)./besselj(m,k0*n_2*radius)+(-1i).^m.*besselj(m,k0*n_1*radius)./besselj(m,k0*n_2*radius);
% 
% sigma =ones(size(m))*2;
% sigma(1)=1;

alphaTE = -1i*8*b_TE(2)/(params.omega*params.mu0);
end
    
    
    
    
    
    