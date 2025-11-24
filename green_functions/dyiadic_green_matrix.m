function [Gxx,Gxy,Gyx,Gyy] = dyiadic_green_matrix(source,test,n,k0,c,OMEGA, params)


xsource = source(1,:);
ysource = source(2,:);
xtest = test(1,:);
ytest = test(2,:);

b = ysource;
a = xsource;
y = ytest;
x = xtest;

OMEGAc = OMEGA/c;


rot = exp(1i*k0*OMEGA/c*(xsource.*ytest-ysource.*xtest));
rr =sqrt((xtest-xsource).*(xtest-xsource)+(ytest-ysource).*(ytest-ysource));
%% first matrix
%dyy
d_yy_1 = 1i/4.*k0.*n.*besselh(1,k0.*n.*rr).*((2.*(y-b).^2./(rr.^3)-1./rr)+1i.*k0.*OMEGAc.*(a+x).*(b-y)./rr).*rot;
d_yy_2 = 1i/4.*(-k0^2.*n^2.*(b-y).^2./(rr.^2).*besselh(0,k0.*n.*rr)+1i.*k0.*OMEGAc.*x.*1i.*k0.*OMEGAc.*a.*besselh(0,k0.*n.*rr)).*rot;

dyy = d_yy_1 + d_yy_2;

%dxx
d_xx_1 = 1i/4*k0*n.*besselh(1,k0*n*rr).*((2*(x-a).^2./(rr.^3)-1./rr)+1i*k0*OMEGAc*(x-a).*(y+b)./rr).*rot;
d_xx_2 = 1i/4*(-k0^2*n^2*(a-x).^2./(rr.^2).*besselh(0,k0*n.*rr)+1i*k0*OMEGAc.*y*1i*k0*OMEGAc.*b.*besselh(0,k0*n.*rr)).*rot;

dxx = d_xx_1+d_xx_2;
%dxb (x,yprime)
d_xb_1 = -1i/4*k0*n.*besselh(1,k0*n.*rr).*(2*(b-y).*(a-x)./(rr.^3)-(b-y)./rr.*(1i*k0*OMEGAc.*b)+1i*k0*OMEGAc*(a-x)./rr.*x).*rot;
d_xb_2 = 1i/4.*besselh(0,k0*n.*rr).*((k0*n)^2*(b-y).*(a-x)./(rr.^2)+(1i*k0*OMEGAc).^2.*x.*b-1i*k0*OMEGAc).*rot;

dxb = d_xb_1+d_xb_2;
%dyb (y,xprime)

d_ya_1 = -1i/4*k0*n.*besselh(1,k0*n*rr).*(2*(b-y).*(a-x)./(rr.^3)+(a-x)./rr.*(1i*k0*OMEGAc.*a)-1i*k0*OMEGAc*(b-y).*y./rr).*rot;
d_ya_2 = 1i/4*besselh(0,k0*n*rr).*((k0*n)^2*(b-y).*(a-x)./(rr.^2)+(1i*k0*OMEGAc).^2.*y.*a+1i*k0*OMEGAc).*rot;

dya = d_ya_1+d_ya_2;
%% static 

Gxx = 1i/(params.omega*(params.e0*params.er_out))*dyy;
Gxy = 1i/(params.omega*(params.e0*params.er_out))*dya;
Gyx = 1i/(params.omega*(params.e0*params.er_out))*dxb;
Gyy = 1i/(params.omega*(params.e0*params.er_out))*dxx;

%% 2nd Matrix

dx = -1i/4*(k0*n*(x-a)./(sqrt((x-a).^2 + (y-b).^2  )).*besselh(1,k0*n*sqrt((x-a).^2 + (y-b).^2))+1i*k0*OMEGAc*b.*besselh(0,k0*n*sqrt((x-a).^2 + (y-b).^2))).*rot;
dy = -1i/4*(k0*n*(y-b)./(sqrt((x-a).^2 + (y-b).^2  )).*besselh(1,k0*n*sqrt((x-a).^2 + (y-b).^2))-1i*k0*OMEGAc*a.*besselh(0,k0*n*sqrt((x-a).^2 + (y-b).^2))).*rot;
da = -1i/4*(k0*n*(a-x)./(sqrt((x-a).^2 + (y-b).^2  )).*besselh(1,k0*n*sqrt((x-a).^2 + (y-b).^2))-1i*k0*OMEGAc*y.*besselh(0,k0*n*sqrt((x-a).^2 + (y-b).^2))).*rot;
db = -1i/4*(k0*n*(b-y)./(sqrt((x-a).^2 + (y-b).^2  )).*besselh(1,k0*n*sqrt((x-a).^2 + (y-b).^2))+1i*k0*OMEGAc*x.*besselh(0,k0*n*sqrt((x-a).^2 + (y-b).^2))).*rot;


B11 = -a.*dy - x.*db;
B12 =  x.*da - b.*dy;
B21 = -a.*dx - y.*db;
B22 =  b.*dx + y.*da;

etha_1 = sqrt(params.mr_out/params.er_out);

Gxx = Gxx + B11 * OMEGAc * etha_1;
Gxy = Gxy + B12 * OMEGAc * etha_1;
Gyx = Gyx + B21 * OMEGAc * etha_1;
Gyy = Gyy + B22 * OMEGAc * etha_1;

%% 3rd Matrix




end