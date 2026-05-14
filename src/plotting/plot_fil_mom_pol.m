function plot_fil_mom_pol(results)
% plot_fil_mom_pol  Compare MoM, Filaments, and Polarizability results across
% a sweep of OMEGA values.
%
% Input:
%   results  : cell array { i, j } of structs returned by eval_TM_results.
%              The j dimension is taken to be the OMEGA sweep.

[ni, nj] = size(results);

% Use the first slot as a representative for plot titles
params_final = results{1,1}.params;

E_MoM_final = [];
E_fil_final = [];
E_pol_final = [];
OMEGA_final = zeros(1, nj);

for i = 1:ni
    for j = 1:nj
        r = results{i,j};
        if isfield(r.mom, 'mean_E')
            E_MoM_final(:, j) = r.mom.mean_E(:);
        end
        E_fil_final(:, j) = r.filaments.mean_E(:);
        E_pol_final(:, j) = full(r.pol.E_sol(:));
        OMEGA_final(j) = r.params.OMEGA_vec(j) / r.params.omega;
    end
end

ttl = sprintf('\\lambda = %g \\mum, Radius = %.4f\\lambda, \\epsilon = %g', ...
              params_final.lambda*1e6, params_final.radius/params_final.lambda, params_final.er_in);

plot_panel(@abs,  OMEGA_final, E_MoM_final, E_fil_final, E_pol_final, ['abs(E) vs \Omega/\omega, ' ttl]);
plot_panel(@real, OMEGA_final, E_MoM_final, E_fil_final, E_pol_final, ['real(E) vs \Omega/\omega, ' ttl]);
plot_panel(@imag, OMEGA_final, E_MoM_final, E_fil_final, E_pol_final, ['imag(E) vs \Omega/\omega, ' ttl]);

% Centered-abs panel (subtract abs at center index of each curve)
figure;
for ll = 1:size(E_MoM_final, 1)
    mid = round(size(E_MoM_final, 2) / 2);
    plot(OMEGA_final, abs(E_MoM_final(ll,:)) - abs(E_MoM_final(ll,mid)) + abs(E_pol_final(ll,mid)), 'x-'); hold on;
    plot(OMEGA_final, abs(E_fil_final(ll,:)) - abs(E_fil_final(ll,mid)) + abs(E_pol_final(ll,mid)), 's-');
    plot(OMEGA_final, abs(E_pol_final(ll,:)), 'd-');
end
legend({'MoM', 'Fil', 'Pol'});
title(['centered abs(E) vs \Omega/\omega, ' ttl]);
xlabel('\Omega/\omega'); grid on; grid minor;

figure;
scatter(params_final.sca_x*1e6, params_final.sca_y*1e6, 500, abs(E_MoM_final(:,1)), '.');
colormap('jet'); grid on; grid minor;
xlabel('\lambda'); ylabel('\lambda');
title('Scatterers Array');

end

function plot_panel(op, x, A, B, C, ttl)
figure;
for ll = 1:size(A, 1)
    plot(x, op(A(ll,:)), 'x-'); hold on;
    plot(x, op(B(ll,:)), 's-');
    plot(x, op(C(ll,:)), 'd-');
end
legend({'MoM', 'Fil', 'Pol'});
title(ttl); xlabel('\Omega/\omega');
grid on; grid minor;
end
