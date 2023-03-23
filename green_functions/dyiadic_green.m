function [Gx,Gy] = dyiadic_green(source,test,params, inside)

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
Iz = params.Iz;
omega = params.omega;

xsource = source(1,:);
ysource = source(2,:);
xtest = test(1,:);
ytest = test(2,:);

rot = exp(1i*k0*OMEGA/c*(xsource.*ytest-ysource.*xtest));

rr =sqrt((xtest-xsource).*(xtest-xsource)+(ytest-ysource).*(ytest-ysource));

%% old version
fact = 1/(1i*omega*mu);

fact_n = 1i*omega*OMEGA/(c*c);

Gst = 1i/4*besselh(0,1,k0*n*rr);
Gst_d = k0*n*1i/4*besselh(1,1,k0*n*rr);

Gx_1 = ((ysource-ytest)./(rr)).*Gst_d;
Gx_2 = -(xtest-xsource).*Gst.*fact_n;

Gy_1 = (-(xsource-xtest)./(rr)).*Gst_d;
Gy_2 = -(ytest-ysource).*Gst.*fact_n;

Gx = Gx_1 + Gx_2;
Gy = Gy_1 + Gy_2;

Gx = Gx .* rot * fact ;
Gy = Gy .* rot * fact ;

end