
function [G] = scalar_green(source,test,n,k0,c,OMEGA,params)

xsource = source(1,:);
ysource = source(2,:);
xtest = test(1,:);
ytest = test(2,:);

rot = exp(1i*k0*OMEGA/c*(xsource.*ytest-ysource.*xtest));
rr = sqrt((xtest-xsource).*(xtest-xsource)+(ytest-ysource).*(ytest-ysource));

Gst = 1i/4*besselh(0,1,k0*n*rr);

G = Gst.*rot;

end