function params = set_bif_param(params, bif_param, param_value)
    % This function updates the selected bifurcation parameter in the params structure
    switch bif_param
        case 'b'
            params.b = param_value;  % Update the 'b' parameter
        case 'r1'
            params.r1 = param_value;  % Update the 'r1' parameter
        case 'r2'
            params.r2 = param_value;  % Update the 'r2' parameter
        case 'f'
            params.f = param_value;  % Update the 'f' parameter
        case 'K'
            params.K = param_value;  % Update the 'K' parameter
        otherwise
            disp('No such bifurcation parameter');
    end
end
