function X = Polynomial_2d(State_Grid, D)
% POLYNOMIAL_2D: Generates the basis matrix X for a 2D ordinary polynomial 
% approximation up to total degree D. (Used for the OLS regression X*v = k')
%
%   Inputs:
%       State_Grid - Matrix of grid points [k, z], size M x 2.
%       D          - Total degree of the polynomial.
%
%   Output:
%       X          - Basis matrix, size M x N_basis.

    M = size(State_Grid, 1);
    k_vec = State_Grid(:, 1); 
    z_vec = State_Grid(:, 2); 
    
    N_basis = (D + 1) * (D + 2) / 2;
    X = zeros(M, N_basis);
    
    col_idx = 1;
    
    for d_total = 0:D 
        for p = 0:d_total 
            q = d_total - p; 
            
            % Basis term: k^p * z^q
            X(:, col_idx) = (k_vec .^ p) .* (z_vec .^ q);
            col_idx = col_idx + 1;
        end
    end

end