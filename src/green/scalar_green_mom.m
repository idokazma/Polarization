function [G] = scalar_green_mom(source,test,n,k0,c,OMEGA,~)
% scalar_green_mom  Variant of scalar_green that takes explicit
% (n, k0, c, OMEGA) instead of the params struct. Used by MoM.m.
% The trailing argument is accepted but ignored for backward compatibility.

xsource = source(1,:);
ysource = source(2,:);
xtest   = test(1,:);
ytest   = test(2,:);

rot = exp(1i*k0*OMEGA/c*(xsource.*ytest - ysource.*xtest));
rr  = sqrt((xtest-xsource).^2 + (ytest-ysource).^2);

Gst = 1i/4*besselh(0,1,k0*n*rr);

G = Gst.*rot;

end
