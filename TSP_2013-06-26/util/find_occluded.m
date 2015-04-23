function [IMG] = find_occluded(IMG)

[IMG.K, IMG.label, IMG.SP, ~, ~, ~, ~, ~] = localonly_move(IMG,0);



flow = reshape([IMG.SP(1:IMG.prev_K).v],[2, IMG.prev_K])*sqrt(IMG.hyper.op_Sigma(1));
u = flow(2,:);
v = flow(1,:);
[utotal, vtotal] = flow_total(IMG.prev_K,IMG.xdim,IMG.ydim,IMG.prev_indices,u,v);
utotal(utotal>10^5) = 0;
vtotal(vtotal>10^5) = 0;
mask = ones(size(IMG.label));
[~, mask] =warpRIMPORT(utotal, vtotal, mask);
IMG.mask_disocclusion = mask<0.5;

IMG.mask_disocclusion = imdilate(~imdilate(~IMG.mask_disocclusion, [0 1 0; 1 1 1; 0 1 0]), [0 1 0; 1 1 1; 0 1 0]) & IMG.boundary_mask;
IMG.mask_disocclusion = imdilate(IMG.mask_disocclusion, [0 1 0; 1 1 1; 0 1 0]);











labelmap = zeros(IMG.K,1,'uint32');
i = 1;
totali = 1;
mask = IMG.label>=0;

while i<=IMG.K
    if (i>numel(IMG.SP))
        break;
    end
    
    if (~IMG.SP(i).old && (isempty(IMG.SP(i).N) || IMG.SP(i).N==0 || isinf(IMG.SP(i).a_mu(1)) || isnan(IMG.SP(i).a_mu(1))))
        IMG.SP(i) = [];
        IMG.K = IMG.K - 1;

    else
        labelmap(totali) = i;
        i = i + 1;
    end
    totali = totali + 1;
end
IMG.label(mask) = labelmap(IMG.label(mask)+1)-1;

app_mean = reshape([IMG.SP.a_mu], [3, IMG.K]);
pos_mean = reshape([IMG.SP.p_mu], [2, IMG.K]);

if (any(isnan(app_mean(:))) || any(isnan(pos_mean(:))))
    disp('NAN ERROR!');
    error(1);
end

mu_a = bsxfun(@times, app_mean', sqrt(IMG.hyper.oa_Sigma))';
mu_p = bsxfun(@times, pos_mean', sqrt(IMG.hyper.op_Sigma))';
mu = cat(1, mu_p, mu_a);
[covariance, ~] = get_gp_covariance(IMG.label, mu, IMG.cov_var_a, IMG.cov_var_p);

theta_p = reshape([IMG.SP(1:IMG.prev_K).p_mu] - [IMG.SP(1:IMG.prev_K).p_theta], [2, IMG.prev_K]);
theta_p = bsxfun(@times, theta_p', sqrt(IMG.hyper.op_Sigma));

new_flow = covariance(:, 1:IMG.prev_K) * (covariance(1:IMG.prev_K, 1:IMG.prev_K) \ theta_p);
[new_indices, ~] = populate_indices(double(IMG.K), IMG.label);
[new_utotal, new_vtotal] = flow_total(IMG.K,IMG.xdim,IMG.ydim,new_indices,new_flow(:,2),new_flow(:,1));
new_utotal(new_utotal>10^5) = 0;
new_vtotal(new_vtotal>10^5) = 0;

mask = ones(size(IMG.label));
[~, mask] =warpRIMPORT(new_utotal, new_vtotal, mask);
IMG.mask_occlusion = mask>1.5;

IMG.mask_occlusion = imdilate(~imdilate(~IMG.mask_occlusion, [0 1 0; 1 1 1; 0 1 0]), [0 1 0; 1 1 1; 0 1 0]) & IMG.boundary_mask;
IMG.D_mask_dilated = imdilate(IMG.mask_occlusion, [0 1 0; 1 1 1; 0 1 0]);

% % fix_stuff
% temp_changed = IMG.SP_changed;
% IMG.SP_changed(:) = false;
% uvals = unique(IMG.label(IMG.D_mask_dilated));
% uvals = uvals(uvals>=0) + 1;
% IMG.SP_changed(uvals) = true;
% [IMG.K, IMG.label, IMG.SP, ~, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy, ~] = localonly_move(IMG,1000);
% IMG.SP_changed(:) = false;
% IMG.SP_changed = IMG.SP_changed | temp_changed;


figure(6);
imagesc(IMG.mask_occlusion);
figure(7);
imagesc(IMG.mask_disocclusion);

a=1