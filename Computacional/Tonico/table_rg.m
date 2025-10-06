function table_rg(x_output, phi_values, rg)

% Displays a table in the console with the values of the scalar field Phi
% at the outlet for different x positions and for different values of the ratio ρ/Γ
%
% Inputs:
%   - x_output   : Vector with x positions at the outlet
%   - phi_values : Matrix with Phi values 
%                  (each column corresponds to a different ρ/Γ)
%   - rg         : Vector with the tested ρ/Γ values
%---------------------------------------------------------------------------

    % Encabezado
    fprintf('\n%-10s', 'x-position');
    for j = 1:length(rg)
        fprintf('%16s', ['ρ/Γ = ', num2str(rg(j), '%.0e')]);
    end
    fprintf('\n');

    % Filas
    for i = 1:length(x_output)
        fprintf('%-10.1f', x_output(i));
        for j = 1:length(rg)
            fprintf('%16.4f', phi_values(i, j));
        end
        fprintf('\n');
    end
end

