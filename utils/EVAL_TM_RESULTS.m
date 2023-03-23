%% EVAL_ALL_RESULTS


original_radius = params.radius;
[E_SOL_st_MoM, mean_E_MoM,mean_I_MoM,effective_radius, alpha_MoM, mean_E0_MoM] = MoM(params);
alpha_MoM
% mean_E_MoM
% mean_I_MoM
close all;
params.radius = effective_radius;
% params.radius = original_radius*0.9874;

[Pvec, E_SOL_st_POL, alpha_Pol] = RotatingArray_2D_TM(params, params.E_inc_z, 0);
alpha_Pol
% E_SOL_st_POL
% Pvec
if (params.is_plane_wave == 1 && length(params.sca_x)== 1 && (params.sca_x)== 0 && (params.sca_y)== 0)
[E_SOL_st_mie, mean_E_mie] = Mie_Series_TM(params,params.E_inc_z);
mean_E_mie = mean(mean_E_mie)
        fprintf('done MIE\n')
end
% params.radius = original_radius * 0.9819;

[E_SOL_st_filaments_mul, mean_E_fil,mean_I_fil, alpha_Fil,field_current_mat] = filaments_TM_multiple(params,0);
alpha_Fil


     
        
        
        
   PARAMETERS_PLANE{i,j,t} = params;
try
    E_SOL_st_mie_PLANE{i,j} = E_SOL_st_mie;
    E_SOL_st_mie_MEAN_PLANE{i,j} = mean_E_mie;
catch
    E_SOL_st_mie_PLANE{i,j} = -1; 
end

E_SOL_st_filaments_mul_PLANE{i,j,t} = E_SOL_st_filaments_mul;
E_SOL_st_filaments_mul_E_PLANE{i,j,t} = mean_E_fil;
E_SOL_st_filaments_mul_I_PLANE{i,j,t} = mean_I_fil;
E_SOL_st_filaments_mul_alpha{i,j,t} = alpha_Fil;
    field_current_mat_output {i,j,t} = field_current_mat;



    
try
E_SOL_st_MoM_PLANE{i,j,t} = E_SOL_st_MoM;
E_SOL_st_MoM_MEAN_PLANE{i,j,t} = mean_E_MoM;
E_SOL_st_MoM_MEAN_PLANE_E0{i,j,t} = mean_E0_MoM;
E_SOL_st_MoM_I_MEAN_PLANE{i,j,t} = mean_I_MoM;
E_SOL_st_alpha_MoM{i,j,t} = alpha_MoM;
catch
E_SOL_st_MoM_PLANE{i,j,t} = -1;
E_SOL_st_MoM_MEAN_PLANE{i,j,t} = -1;
E_SOL_st_MoM_MEAN_PLANE_E0{i,j,t} = -1;
E_SOL_st_MoM_I_MEAN_PLANE{i,j,t} = -1;
E_SOL_st_alpha_MoM{i,j,t} = -1;    
end


E_SOL_st_POL_PLANE{i,j,t} = E_SOL_st_POL;
E_SOL_st_POL_I_PLANE{i,j,t} = Pvec;
E_SOL_st_alpha_POL{i,j,t} = alpha_Pol;

