
function [G] = scalar_green(source,test,params, inside)

if (inside)
n = params.n_in;
mu = (params.mu0*params.mr_in);
else
n = params.n_out;
mu = (params.mu0*params.mr_out); 
end

c = params.c;
OMEGA = params.OMEGA;
k0 = params.k0;
Iz = 1;
omega = params.omega;

xsource = source(1,:);
ysource = source(2,:);
xtest = test(1,:);
ytest = test(2,:);

rot = exp(1i*k0*OMEGA/c*(xsource.*ytest-ysource.*xtest));
rr =sqrt((xtest-xsource).*(xtest-xsource)+(ytest-ysource).*(ytest-ysource));

fact = -omega*mu*Iz/4;
Gst = fact* besselh(0,1,k0*n*rr);

G = Gst.*rot;

end