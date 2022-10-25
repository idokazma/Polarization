%% Mie Series TM Calculator
function [full_fields,mean_E_mie] = Mie_Series_TM(params,E_inc_z)

[theta,rho] = cart2pol(params.X,params.Y);
x = params.x;
y = params.y;


E_0 = E_inc_z(params.X==0 & params.Y==0);

n_1 = params.n_out;
n_2 = params.n_in;

k0 = params.k0;
radius = params.radius;

etha_1 = sqrt(params.mr_out/params.er_out);
etha_2 = sqrt(params.mr_in/params.er_in);

m=0:params.max_m;
tic
%% b_tm calc
b_m_tm_up_A = etha_2*(1/2*(besselj(m-1,k0*n_1*radius)-besselj(m+1,k0*n_1*radius))).*besselj(m,k0*n_2*radius);
b_m_tm_up_B = etha_1*(1/2*(besselj(m-1,k0*n_2*radius)-besselj(m+1,k0*n_2*radius))).*besselj(m,k0*n_1*radius);

b_m_tm_down_C = etha_2*(1/2*(besselh(m-1,1,k0*n_1*radius)-besselh(m+1,1,k0*n_1*radius))).*besselj(m,k0*n_2*radius);
b_m_tm_down_D = etha_1*(1/2*(besselj(m-1,k0*n_2*radius)-besselj(m+1,k0*n_2*radius))).*besselh(m,1,k0*n_1*radius);

b_TM = ((-1)*(-1i).^m).*(b_m_tm_up_A-b_m_tm_up_B)./(b_m_tm_down_C-b_m_tm_down_D);

%% c_tm calc

c_TM = b_TM.*besselh(m,1,k0*n_1*radius)./besselj(m,k0*n_2*radius)+(-1i).^m.*besselj(m,k0*n_1*radius)./besselj(m,k0*n_2*radius);

sigma =ones(size(m))*2;
sigma(1)=1;

E_mid=zeros(size(params.X));
E_sca = zeros(size(params.X));
E_final = E_mid;

H_rho_mid = zeros(size(params.X));
H_theta_mid = zeros(size(params.X));

H_x_mid = zeros(size(params.X));
H_y_mid = zeros(size(params.X));

for(i=1:length(m))
E_mid(rho>radius) =  E_mid(rho>radius)+sigma(i)*b_TM((i)).*besselh(m(i),1,k0*n_1.*rho(rho>radius)).*cos(m(i)*theta(rho>radius));
 E_mid(rho<=radius) =  E_mid(rho<=radius)+sigma(i)*c_TM(i).*besselj(m(i),k0*n_2.*rho(rho<=radius)).*cos(m(i)*theta(rho<=radius));

  H_rho_mid(rho<=radius) =  H_rho_mid(rho<=radius)+sigma(i)*c_TM(i).*besselj(m(i),k0*n_2.*rho(rho<=radius)).*(-m(i)*sin((m(i)*theta(rho<=radius))))./(rho(rho<=radius));
  H_theta_mid(rho<=radius) =  H_theta_mid(rho<=radius)+sigma(i)*c_TM(i).*(0.5*k0*n_2*(besselj(m(i)-1,k0*n_2.*rho(rho<=radius))-besselj(m(i)+1,k0*n_2.*rho(rho<=radius))))...
     .*cos(m(i)*theta(rho<=radius));

 H_rho_mid(rho>radius) =  H_rho_mid(rho>radius)+sigma(i)*b_TM((i)).*besselh(m(i),1,k0*n_1.*rho(rho>radius)).*(-m(i)*sin((m(i)*theta(rho>radius))))./(rho(rho>radius));
 H_theta_mid(rho>radius) =  H_theta_mid(rho>radius)+sigma(i)*b_TM((i)).*(0.5*k0*n_1*(besselh(m(i)-1,1,k0*n_1.*rho(rho>radius))-besselh(m(i)+1,1,k0*n_1.*rho(rho>radius)))).*cos(m(i)*theta(rho>radius));

 
end
H_rho_mid(rho<=radius) = H_rho_mid(rho<=radius)./(1i*params.omega*params.mr_in*params.mu0).*E_0;
H_rho_mid(rho>radius) = H_rho_mid(rho>radius)./(1i*params.omega*params.mr_out*params.mu0).*E_0;

H_theta_mid(rho<=radius) = H_theta_mid(rho<=radius)./(1i*params.omega*params.mr_in*params.mu0).*E_0;
H_theta_mid(rho>radius) = H_theta_mid(rho>radius)./(1i*params.omega*params.mr_out*params.mu0).*E_0;


H_theta_mid = -H_theta_mid;

H_x_mid = cos(theta).*H_rho_mid -sin(theta).*H_theta_mid;
H_y_mid = sin(theta).*H_rho_mid +cos(theta).*H_theta_mid;

H_x_final = H_x_mid;
H_x_final(rho>radius) = H_x_final(rho>radius) + (params.H_inc_x(rho>radius));

H_y_final = H_y_mid;
H_y_final(rho>radius) = H_y_final(rho>radius) + (params.H_inc_y(rho>radius));

E_final(rho>radius) = E_0*E_mid(rho>radius)+(E_inc_z(rho>radius));
E_final(rho<=radius)= E_0*E_mid(rho<=radius);

full_fields = {E_final,H_x_final,H_y_final};

mean_E_mie = mean(E_final(rho<=radius));
toc

E_sca(rho>radius) = E_0*E_mid(rho>radius);
E_sca(rho<=radius)= E_0*E_mid(rho<=radius) - (E_inc_z(rho<=radius));

% figure; imagesc(x,y,abs(E_mid)); set(gca,'YDir','normal'); axis equal; title ('ENG Mie Series');


E_SOL_st_mie = E_final;

end