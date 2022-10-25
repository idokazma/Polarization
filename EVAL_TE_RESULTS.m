%% EVAL_ALL_RESULTS



if (params.is_plane_wave == 1 && length(params.sca_x)== 1 && (params.sca_x)== 0 && (params.sca_y)== 0)
% [H_SOL_st_mie, H_sca_mie, Ex_mie, Ey_mie, Ex_sca_mie, Ey_sca_mie, alpha_TE_mie] = Mie_Series_TE(params,params.H_inc_z);
% alpha_TE_mie

% mean_H_mie = mean(mean_H_mie)
        fprintf('done MIE\n')
end
[H_SOL_st_filaments_mul, H_SOL_st_filaments_sca, Ex_sol_fil, Ey_sol_fil,Ex_sol_final_sca,Ey_sol_final_sca, alpha_mat_fil,field_current_mat] = filaments_TE_multiple(params,0);
alpha_mat_fil        
fprintf ('done FIL\n')


     
        
        
        
  PARAMETERS_PLANE{i,j,t} = params;
  alpha_te_fil_output {i,j,t} = alpha_mat_fil;
    field_current_mat_output {i,j,t} = field_current_mat;

% try
%     H_SOL_st_mie_PLANE{i,j} = H_SOL_st_mie;
% %     E_SOL_st_mie_MEAN_PLANE{i,j} = mean_E_mie;
% catch
%     H_SOL_st_mie_PLANE{i,j} = -1; 
% end
% 
% E_SOL_st_filaments_mul_PLANE{i,j} = E_SOL_st_filaments_mul;
% E_SOL_st_filaments_mul_I_PLANE{i,j} = mean_E_fil;
% 
% E_SOL_st_MoM_PLANE{i,j} = E_SOL_st_MoM;
% E_SOL_st_MoM_MEAN_PLANE{i,j} = mean_E_MoM;
% E_SOL_st_MoM_I_MEAN_PLANE{i,j} = mean_I_MoM;
% 
% E_SOL_st_POL_PLANE{i,j} = E_SOL_st_POL;
% E_SOL_st_POL_I_PLANE{i,j} = Pvec;
toc

fprintf ('~~~~ DONE eps %d, radius %d ~~~~ \n', i/length(eps_vec)*100, j/length(radius_vec)*100)