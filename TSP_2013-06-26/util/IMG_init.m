function IMG = IMG_init(img, params)

IMG.cov_var_a = params.cov_var_a;
IMG.cov_var_p = params.cov_var_p;
%IMG.prev_iid_var = params.prev_iid_var;
IMG.alive_dead_changed = true;

% Build IMG structure

% 1. image statistics
IMG.oxdim = size(img,1);
IMG.oydim = size(img,2);
%IMG.label = zeros([IMG.xdim,IMG.ydim]);

% % RGB
% IMG.data(3,:) = reshape(img(:,:,1),1,num_pix);
% IMG.data(4,:) = reshape(img(:,:,2),1,num_pix);
% IMG.data(5,:) = reshape(img(:,:,3),1,num_pix);

N = IMG.oxdim*IMG.oydim;
IMG.area = N/params.K;
IMG.area_var = params.area_var;


IMG.w = 2*round(2*sqrt(IMG.area/pi));
IMG.w = 2*round(2*sqrt(IMG.area/pi));

% IMG.w = 0;

IMG.xdim = IMG.oxdim + 2*IMG.w;
IMG.ydim = IMG.oydim + 2*IMG.w;

IMG.boundary_mask = false(IMG.xdim, IMG.ydim);
IMG.boundary_mask(IMG.w+1:end-IMG.w, IMG.w+1:end-IMG.w) = true;


% IMG.hyper.p_Sigma = repmat(area/3, [1,2]);
% IMG.hyper.p_Sigma = [50 50];
% IMG.hyper.p_Delta = [100 100];
% IMG.hyper.a_Sigma = [100 10 10]*25;
% IMG.hyper.a_Delta = [1000 100 100]/100;


IMG.log_alpha = params.alpha * IMG.area;
IMG.log_beta = params.beta * IMG.area;
Sigma = IMG.area^2 / (-12.5123*IMG.log_alpha);


IMG.hyper.p_Sigma = [Sigma Sigma];
IMG.hyper.p_Delta = [Sigma*2 Sigma*2]*params.deltap_scale;
% IMG.hyper.a_Sigma = [Sigma*2 Sigma Sigma]*15;
IMG.hyper.a_Sigma = [Sigma*2 Sigma Sigma]*params.K/100;
IMG.hyper.a_Delta = [Sigma*20 Sigma*10 Sigma*10]/params.deltaa_scale;

r = sqrt(IMG.area / pi);
IMG.dummy_log_prob = (-0.5 * r^2/IMG.hyper.p_Sigma(1)) - log(2*pi .* IMG.hyper.p_Sigma(1));
% IMG.dummy_log_prob = -1000;
% IMG.dummy_log_prob = 10000;
% IMG.dummy_log_prob = IMG.dummy_log_prob*1.5;

IMG.hyper.p_theta = [0 0];
IMG.hyper.a_theta = [0 0 0];
IMG.hyper.op_Sigma = IMG.hyper.p_Sigma;
IMG.hyper.oa_Sigma = IMG.hyper.a_Sigma;

IMG.data = SP_img2data(img, IMG.w);
IMG.data = rescale_data(IMG.data, IMG.hyper.op_Sigma, IMG.hyper.oa_Sigma);

IMG.hyper.p_theta = IMG.hyper.p_theta ./ sqrt(IMG.hyper.p_Sigma);
IMG.hyper.p_Delta = IMG.hyper.p_Delta ./ (IMG.hyper.p_Sigma);
IMG.hyper.p_Sigma(:) = 1;

IMG.hyper.a_theta = IMG.hyper.a_theta ./ sqrt(IMG.hyper.a_Sigma);
IMG.hyper.a_Delta = IMG.hyper.a_Delta ./ (IMG.hyper.a_Sigma);
IMG.hyper.a_Sigma(:) = 1;


% -1 1
% logN prior
%IMG.alpha = -1.3;
%IMG.epsilon = 1;

% DP-like prior
%IMG.alpha = 1;

%IMG.alpha = exp(-(N/K-sqrt(K)-2.6));
%IMG.alpha = exp(-(N/K - 80));
%IMG.alpha = -N/K * 0.5;
%IMG.alpha = -100;
% should be x+30 ~ N/K
% x = N/K-30

% IMG.log_alpha = -log((N/K + 870) / 844) *715;
% IMG.log_alpha = -(N / (K-135) * 0.33 - 17);
% 
% IMG.log_alpha = IMG.log_alpha + 36;
% 
% IMG.log_beta = IMG.log_alpha + 600;
% IMG.log_beta = IMG.log_alpha + 50;

% IMG.log_alpha = -30;
% IMG.log_alpha = -alpha;
% IMG.log_beta = -50;
% IMG.log_beta = -75;

%IMG.alpha = 0.5;
%IMG.alpha = 1;
%disp(['alpha=' num2str(IMG.alpha)]);
%IMG.alpha = 0;

% 2. SuperPixel statistics
[IMG.label, IMG.K] = random_init(IMG.xdim,IMG.ydim,IMG.w,round(params.K*params.Kpercent));
IMG.label(~IMG.boundary_mask) = 0; % this will be -1 in the end

% every pixel has it's own cluster
% IMG.label(:) = 1:N;
% gridded
% IMG.label(:) = floor(yg/10)*4 + floor(xg/10) + 1;
% IMG.label(9:13,9:13) = max(IMG.label(:))+1;
IMG.K = double(max(IMG.label(:)));

%IMG.label(:) = 1:numel(IMG.label(:));
%IMG.K = numel(IMG.label(:));
IMG.label = IMG.label-1;


%IMG.SP = cell(1,N);
IMG.SP = struct(...
    'p_Delta',{}, 'p_theta',{}, 'p_mu',{}, 'p_Sigma',{}, ...
    'a_Delta',{}, 'a_theta',{}, 'a_mu',{}, 'a_Sigma',{}, ...
    'N',{}, 'UID',{}, 'v',{}, 'prev_v', {}, 'old', {});


% for i = 1: IMG.K
%     %position
%     IMG.SP(i).p_kappa = IMG.hyper.p_kappa;
%     IMG.SP(i).p_nu = IMG.hyper.p_nu;
%     IMG.SP(i).p_Delta = IMG.hyper.p_Delta;
%     IMG.SP(i).p_theta = IMG.hyper.p_theta;
%     IMG.SP(i).p_mu = zeros(1,2);
%     IMG.SP(i).p_Sigma = zeros(2,2);
%     %appearance
%     IMG.SP(i).a_kappa = IMG.hyper.a_kappa;
%     IMG.SP(i).a_nu = IMG.hyper.a_nu;
%     IMG.SP(i).a_Delta = IMG.hyper.a_Delta;
%     IMG.SP(i).a_theta = IMG.hyper.a_theta;
%     IMG.SP(i).a_mu = zeros(1,3);
%     IMG.SP(i).a_Sigma = zeros(3,3);
% 
%     %2.2 Pixel statistics
%     IMG.SP(i).N = 0;
% end



% 3. Topology Table Look Up
load topology_tables;
IMG.T4Table = tc.T4Table;

IMG.max_UID = zeros(1,1,'uint64');

IMG.SP_changed = true(1, IMG.xdim*IMG.ydim);
