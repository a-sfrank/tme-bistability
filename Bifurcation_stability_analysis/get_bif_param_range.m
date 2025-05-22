% This function returns the range of values for the bifurcation parameter
function param_range = get_bif_param_range(bif_param)
    switch bif_param
        case 'b'
            param_range = 1.0:0.01:4.5;  % Example for 'b' parameter range
        case 'r1'
            param_range = 0.5:0.1:1.5;  % Example for 'r1' parameter range
        case 'r2'
            param_range = 0.5:0.1:1.5;  % Example for 'r2' parameter range
        case 'f'
            param_range = 0.1:0.2:2.4;  % Example for 'f' parameter range
        case 'K'
            param_range = 0.4:0.1:1.1;  % Example for 'K' parameter range
        otherwise
            disp('No such bifurcation parameter');
            param_range = [];  % Return an empty array if parameter is invalid
    end
end
