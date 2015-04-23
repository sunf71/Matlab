function IMG = IMG_prop(img,vx,vy,IMG, params)


% delete the SPs that are only int he boundary
unique_vals_image = unique(IMG.label(IMG.boundary_mask));
unique_vals_boundary = unique(IMG.label(~IMG.boundary_mask));
boundarySPs = setdiff(unique_vals_boundary, unique_vals_image);
IMG.label(ismember(IMG.label, boundarySPs)) = -1;
[IMG.K, IMG.label, IMG.SP, ~, ~, ~, ~, ~] = localonly_move(IMG,0);


% Build IMG structure

% IMG.lambda = 10;
% IMG.lambda_sigma = 25;

% 1. image statistics
IMG.data = SP_img2data(img, IMG.w);
IMG.data = rescale_data(IMG.data, IMG.hyper.op_Sigma, IMG.hyper.oa_Sigma);



%IMG.alpha = 0.5;
%disp(['alpha=' num2str(IMG.alpha)]);
%IMG.alpha = 0;

%
% DW: hard for later tracking...
% delete unused super pixels and change the label as well
labelmap = zeros(IMG.K,1,'uint32');
i = 1;
totali = 1;
mask = IMG.label>=0;

while i<=IMG.K
    if (isempty(IMG.SP(i).N) || IMG.SP(i).N==0)
        IMG.SP(i) = [];
        IMG.K = IMG.K - 1;
    else
        labelmap(totali) = i;
        i = i + 1;
    end
    totali = totali + 1;
end
IMG.label(mask) = labelmap(IMG.label(mask)+1)-1;
%}

meanp = reshape([IMG.SP(:).p_mu], [2,numel(IMG.SP)]);
meanx = meanp(1,:);
meany = meanp(2,:);

%label = SP_prop_init(IMG.K,IMG.label,meanx,meany,sp_vx+randn(1,IMG.K),sp_vy+randn(1,IMG.K));

IMG.prev_label = IMG.label;
IMG.prev_K = IMG.K;
IMG.prev_app_mean = reshape([IMG.SP.a_mu], [3, IMG.prev_K]);
IMG.prev_pos_mean = reshape([IMG.SP.p_mu], [2, IMG.prev_K]);



% get the previous precision
mu_a = bsxfun(@times, IMG.prev_app_mean', sqrt(IMG.hyper.oa_Sigma))';
mu_p = bsxfun(@times, IMG.prev_pos_mean', sqrt(IMG.hyper.op_Sigma))';
mu = cat(1, mu_p, mu_a);
[IMG.prev_covariance, IMG.prev_precision] = get_gp_covariance(IMG.label, mu, IMG.cov_var_a, IMG.cov_var_p, IMG.hyper.p_Delta(1));
% IMG.prev_covariance = 100000*eye(size(IMG.prev_covariance,1));
% IMG.prev_precision = eye(size(IMG.prev_covariance,1))/100000;

vx_extended = holdpad(vx, size(IMG.label,1),size(IMG.label,2));
vy_extended = holdpad(vy, size(IMG.label,1),size(IMG.label,2));

[IMG.prev_indices, ~] = populate_indices(double(IMG.prev_K), IMG.prev_label);

sp_v = zeros(2, numel(IMG.SP));
sp_x = zeros(1, numel(IMG.SP));
sp_y = zeros(1, numel(IMG.SP));
for i = 1: IMG.K
%     IMG.SP(i).p_theta = IMG.SP(i).p_mu;
    indices_k = IMG.prev_indices(i).all;
    IMG.SP(i).a_theta = IMG.SP(i).a_mu;
    
    vxi = mean(vx_extended(indices_k));
    vyi = mean(vy_extended(indices_k));

    x = round((IMG.SP(i).p_mu(1))*sqrt(IMG.hyper.op_Sigma(1)))+1;
    y = round((IMG.SP(i).p_mu(2))*sqrt(IMG.hyper.op_Sigma(2)))+1;
    xi = max(min(x,IMG.xdim-2*IMG.w),1);
    yi = max(min(y,IMG.ydim-2*IMG.w),1);
%     IMG.SP(i).v(1:2) = 0;
    
    %2.2 Pixel statistics
    IMG.SP(i).N = 0;
    IMG.SP(i).old = true;
    
    IMG.SP(i).p_theta = IMG.SP(i).p_mu;% + [vx(xi,yi), vy(xi,yi)] ./ sqrt(IMG.hyper.op_Sigma);
    
    if (isempty(IMG.SP(i).prev_v) || all(IMG.SP(i).prev_v==0))
        IMG.SP(i).prev_v(:) = [vxi, vyi] ./ sqrt(IMG.hyper.op_Sigma);
    else
        IMG.SP(i).prev_v(:) = IMG.SP(i).v(:);
    end
    
%     vxi = IMG.SP(i).prev_v(1) * sqrt(IMG.hyper.op_Sigma(1));
%     vyi = IMG.SP(i).prev_v(2) * sqrt(IMG.hyper.op_Sigma(2));
    
    IMG.SP(i).v = [vx(xi,yi), vy(xi,yi)] ./ sqrt(IMG.hyper.op_Sigma);
    
    sp_v(1,i) = vxi;
    sp_v(2,i) = vyi;
    
    sp_x(i) = vxi + x;
    sp_y(i) = vyi + y;
end
meanx = meanx *sqrt(IMG.hyper.op_Sigma(1));
meany = meany *sqrt(IMG.hyper.op_Sigma(2));


%sp_v = reshape([IMG.SP(:).v], [2,numel(IMG.SP)]);
sp_vx = sp_v(1,:);
sp_vy = sp_v(2,:);


% sp_x = max(min(round(sp_x),IMG.xdim),1);
% sp_y = max(min(round(sp_y),IMG.ydim),1);
% centers = sub2ind([IMG.xdim, IMG.ydim], sp_x, sp_y);
% label = uint32(random_initIMPORT(IMG.xdim,IMG.ydim, centers))-1;


%prev_cov = cellfun(@inv,{IMG.SP.p_Sigma},'UniformOutput',false);
%IMG.prev_pos_prec = reshape([prev_cov{1:end}], [2, 2, IMG.prev_K]);
label = SP_prop_init(IMG.K,IMG.label,meanx,meany,sp_vx,sp_vy,IMG.boundary_mask);
% label(label>=IMG.prev_K & ~IMG.boundary_mask) = -1;
% label = Util_labelconnect(label);
% 
% %sfigure(3);imagesc(label)
% mask = label==IMG.K;
% cc = bwlabel(mask,4);
% label(mask) = label(mask) + uint32(cc(mask) - 1);
IMG.label = label;
IMG.K = double(max(label(:)))+1;

%[label, K] = random_prop(IMG);
%[IMG.label, IMG.K] = random_prop(IMG);
%IMG.label = IMG.label-1;

IMG.SP_changed(:) = true;




% [IMG] = find_occluded(IMG);