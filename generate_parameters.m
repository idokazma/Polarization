

toplot = 1;

start_time = datestr(now);
disp (start_time);
%declare GRID to calculate the fields in;



params.y = linspace(0,params.len,params.len_n)- params.len/2;
params.x = linspace(0,params.wid,params.wid_n)- params.wid/2;

[params.X,params.Y] = meshgrid(params.x,params.y);

% params.len_n_mom = 201;
% params.wid_n_mom = 201;
% params.len_mom = 2.1*1e-6;
% params.wid_mom = 2.1*1e-6;

params.y_mom = linspace(0,params.len,params.len_n*1)- params.len/2;
params.x_mom = linspace(0,params.wid,params.wid_n*1)- params.wid/2;

%% State the physical parameters fo the problem (epsilon(r), mu(r), radius...)
lambda = params.lambda;

% params.radius = 0.01*1e-6; % [meters]



params.e0=8.8541878128e-12 ; 8.85418781762039*1e-12;
e0 = params.e0;

params.mu0 = 1.25663706212e-6;4*pi*1e-7;

mu0 = params.mu0;

params.mr_in = 1;
mr_in = params.mr_in;
params.mr_out = 1;
mr_out = params.mr_out;
params.er_in = 11.4;
er_in = params.er_in;
params.er_out = 1;
er_out = params.er_out;

params.c=1/sqrt(params.e0*params.mu0);
c =params.c;

params.f = params.c/params.lambda;
f = params.f;

params.k0 = 2*pi/params.lambda;
k0 = params.k0;

params.omega = params.f*2*pi;
omega = params.omega;

params.n_in = sqrt(params.mr_in*params.er_in);
params.n_out = sqrt(params.mr_out*params.er_out);



    %    Parameters of Filamets;
    params.R_out=1.2; %params.radius*param.R_out
    params.R_in=0.8;
    
    params.N_filaments=60; %on each side
    params.N_testpoints=200;%must be >=param.N_filaments

    params.max_m = 50;