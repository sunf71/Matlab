function [covariance, precision] = get_gp_covariance(z, mu, cov_var_a, cov_var_p, iid_var)

if (~exist('iid_var', 'var') || isempty(iid_var))
    iid_var = 0;
end

all_unique_z = unique(z(z>=0));
Nz = numel(all_unique_z);
covariance = zeros(Nz);
    
for zi=1:Nz
    covariance(zi,:) = exp(-sum(bsxfun(@minus, mu(1:2,:), mu(1:2,zi)).^2 / (2*cov_var_p))) .* ...
                       exp(-sum(bsxfun(@minus, mu(3:5,:), mu(3:5,zi)).^2 / (2*cov_var_a)));
end

% [V,D] = eig(covariance);
% D(D<0) = 0;
% covariance = real(V)*real(D)*real(V)^-1;
% covariance = 0.5*(covariance + covariance');
precision = (covariance + eye(Nz)*iid_var)^-1;