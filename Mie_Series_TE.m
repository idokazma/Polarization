%% Mie Series TE Calculator
function [H_SOL_st_mie,H_sca,E_x_hat_mie,E_y_hat_mie,E_x_hat_mie_sca,E_y_hat_mie_sca, alpha_TE] = Mie_Series_TE(params,H_inc_z)

[theta,rho] = cart2pol(params.X,params.Y);
x = params.x;
y = params.y;


H_0 = H_inc_z(params.X==0 & params.Y==0);



n_1 = params.n_out;
n_2 = params.n_in;

k0 = params.k0;
radius = params.radius;

etha_1 = sqrt(params.mr_out/params.er_out);
etha_2 = sqrt(params.mr_in/params.er_in);

m=0:params.max_m;
tic
%% b_te calc
b_m_te_up_A = etha_1*(1/2*(besselj(m-1,k0*n_1*radius)-besselj(m+1,k0*n_1*radius))).*besselj(m,k0*n_2*radius);
b_m_te_up_B = etha_2*(1/2*(besselj(m-1,k0*n_2*radius)-besselj(m+1,k0*n_2*radius))).*besselj(m,k0*n_1*radius);

b_m_te_down_C = etha_1*(1/2*(besselh(m-1,1,k0*n_1*radius)-besselh(m+1,1,k0*n_1*radius))).*besselj(m,k0*n_2*radius);
b_m_te_down_D = etha_2*(1/2*(besselj(m-1,k0*n_2*radius)-besselj(m+1,k0*n_2*radius))).*besselh(m,1,k0*n_1*radius);

b_TE = ((-1)*(-1i).^m).*(b_m_te_up_A-b_m_te_up_B)./(b_m_te_down_C-b_m_te_down_D);

%% c_tm calc

c_TE = b_TE.*besselh(m,1,k0*n_1*radius)./besselj(m,k0*n_2*radius)+(-1i).^m.*besselj(m,k0*n_1*radius)./besselj(m,k0*n_2*radius);

sigma =ones(size(m))*2;
sigma(1)=1;

alpha_TE = -1i*8*b_TE(2)/(params.omega*params.mu0)


H_mid=zeros(size(params.X));
H_sca = zeros(size(params.X));
H_final = H_mid;
E_rho_hat = H_mid;
E_theta_hat = H_mid;

OUT_loc = rho>radius;
IN_loc = rho<=radius;

for(i=1:length(m))

 H_mid(OUT_loc) =  H_mid(OUT_loc)+sigma(i)*b_TE((i)).*besselh(m(i),1,k0*n_1.*rho(OUT_loc)).*cos(m(i)*theta(OUT_loc));
 H_mid(IN_loc) =  H_mid(IN_loc)+sigma(i)*c_TE(i).*besselj(m(i),k0*n_2.*rho(IN_loc)).*cos(m(i)*theta(IN_loc));

 E_theta_hat(OUT_loc) = E_theta_hat(OUT_loc) - sigma(i)*b_TE((i)).*(0.5*(besselh(m(i)-1,1,k0*n_1.*rho(OUT_loc))-besselh(m(i)+1,1,k0*n_1.*rho(OUT_loc)))).*cos(m(i).*theta(OUT_loc))*k0*n_1;
 E_rho_hat(OUT_loc) = E_rho_hat(OUT_loc) + sigma(i)*b_TE((i)).*besselh(m(i),1,k0*n_1.*rho(OUT_loc)).*(-1).*m(i).*sin(m(i).*theta(OUT_loc))./rho(OUT_loc);


 E_theta_hat(IN_loc) = E_theta_hat(IN_loc) - sigma(i)*c_TE((i)).*(0.5*(besselj(m(i)-1,k0*n_2.*rho(IN_loc))-besselj(m(i)+1,k0*n_2.*rho(IN_loc)))).*cos(m(i).*theta(IN_loc))*k0*n_2;
 E_rho_hat(IN_loc) = E_rho_hat(IN_loc) + sigma(i)*c_TE((i)).*besselj(m(i),k0*n_2.*rho(IN_loc)).*(-1).*m(i).*sin(m(i).*theta(IN_loc))./rho(IN_loc);

end

E_rho_hat(IN_loc) = H_0*E_rho_hat(IN_loc)*1i/(params.omega*params.e0*params.er_in);
E_theta_hat(IN_loc) = H_0*E_theta_hat(IN_loc)*1i/(params.omega*params.e0*params.er_in);

E_rho_hat(OUT_loc) = H_0*E_rho_hat(OUT_loc)*1i/(params.omega*params.e0*params.er_out);
E_theta_hat(OUT_loc) = H_0*E_theta_hat(OUT_loc)*1i/(params.omega*params.e0*params.er_out);

E_x_hat_mie = cos(theta).*E_rho_hat - sin(theta).*E_theta_hat;
E_y_hat_mie = sin(theta).*E_rho_hat + cos(theta).*E_theta_hat;

E_x_hat_mie(rho>radius) = E_x_hat_mie(rho>radius)+params.E_inc_x(rho>radius);
E_y_hat_mie(rho>radius) = E_y_hat_mie(rho>radius)+params.E_inc_y(rho>radius);

E_x_hat_mie_sca = E_x_hat_mie - params.E_inc_x; 
E_y_hat_mie_sca = E_y_hat_mie - params.E_inc_y; 

H_final(rho>radius) = H_0*H_mid(rho>radius)+(H_inc_z(rho>radius));
H_final(rho<=radius)= H_0*H_mid(rho<=radius);
toc
H_sca = H_final - H_inc_z;
% figure; imagesc(x,y,abs(H_mid)); set(gca,'YDir','normal'); axis equal; title ('ENG Mie Series');

H_SOL_st_mie = H_final;
end