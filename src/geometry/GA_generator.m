function VogelArrayXY = GA_generator(N)
% GA_generator  Generate a Vogel spiral (golden-angle) array of N points in 2D.
% Returns a 2xN array of (x,y) coordinates.

if nargin < 1, N = 100; end

GoldenRatioS = 1.5 + 0.5*sqrt(5);   % square of the golden ratio
alphaGA      = 2*pi/GoldenRatioS;   % the golden angle [rad]
min_dist     = 10;
a            = (1/1.6) * min_dist;

VogelArrayXY = zeros(2, N);
for ip = 1:N
    rp     = sqrt(ip) * a;
    alphap = alphaGA * ip;
    VogelArrayXY(1, ip) = rp*cos(alphap);
    VogelArrayXY(2, ip) = rp*sin(alphap);
end

end
