function [sp_labels] = TSP(K, root, files, dispOn, frames)
%TSP Temporal Superpixel Segmentation.
%   SP_LABELS = TSP(K, ROOT, FILES) returns the label matrix in time and
%   space for the video volume in UINT32. K is the (approximate) number of
%   superpixels per frame. ROOT is the directory to the frames. FILES is a
%   list of the frame images, typically obtained using
%      FILES = dir([ROOT '*.jpg']);
%
%   SP_LABELS = TSP(K, ROOT, FILES, DISPON) supplies an additional flag to
%   display the progress of the algorithm while processing. If omitted or
%   empty, DISPON defaults to true.
%
%   SP_LABELS = TSP(K, ROOT, FILES, DISPON, FRAMES) supplies an additional
%   variable that indicates which frames to process. FRAMES should be in
%   the format of STARTFRAME:ENDFRAME. If omitted or empty, FRAMES defaults
%   to 1:NUMFRAMES.

%   Notes: This version of the code does not reestimate the flow between
%   frames. As noted in the paper, the flow estimation does not do much. If
%   you desire to reestimate the flow, set params.reestimateFlow to be
%   true.
%
%   All work using this code should cite:
%   J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
%      Temporal Superpixels. CVPR 2013.
%
%   Written by Jason Chang and Donglai Wei 2013/06/20

% add the necessary paths
addpath('gui/');
addpath('mex/');
addpath('util/');


params.cov_var_p = 1000;
params.cov_var_a = 100;
params.area_var = 400;
params.alpha = -15;
params.beta = -10;
params.deltap_scale = 1e-3;
params.deltaa_scale = 100;
params.K = K;
params.Kpercent = 0.8;
params.reestimateFlow = false;

if (~exist('dispOn','var') || isempty(dispOn))
    dispOn = true;
end

root_flows = fullfile(root,'TSP_flows/');
if (~exist(root_flows,'dir'))
    mkdir(root_flows);
end

disp('Precomputing all the optical flows...');
for f=2:numel(files)
    im1 = imread(fullfile(root,files(f-1).name));
    im2 = imread(fullfile(root,files(f).name));
    outname = fullfile(root_flows,[files(f).name(1:end-4) '_flow.mat']);
    disp([' -> ' outname]);
    compute_of(im1,im2,outname);
end
disp(' -> Optical flow calculations done');    


flow_files = dir([root_flows '*_flow.mat']);

if (~exist('frames','var') || isempty(frames))
    frames = 1:numel(files);
else
    frames(frames>numel(files)) = [];
end
oim = imread([root files(1).name]);
sp_labels = zeros(size(oim,1), size(oim,2), numel(frames), 'uint32');
frame_it = 0;

disp('Staring Segmentation');
for f=frames
    disp([' -> Frame '  num2str(f) ' / ' num2str(numel(frames))]);

    frame_it = frame_it + 1;
    oim1 = imread([root files(f).name]);

    if (frame_it==1)
        IMG = IMG_init(oim1, params);
    else
        % optical flow returns actual x and y flow... flip it
        vx = zeros(size(oim,1), size(oim,2));
        vy = zeros(size(oim,1), size(oim,2));
        % load the optical flow
        load([root_flows flow_files(f-1).name]);

        vx = -flow.bvx;
        vy = -flow.bvy;
        IMG = IMG_prop(oim1,vy,vx,IMG);
    end

    oim = oim1;

    E = [];
    it = 0;
    IMG.alive_dead_changed = true;
    IMG.SxySyy = [];
    IMG.Sxy = [];
    IMG.Syy = [];
    converged = false;
    while (~converged && it<5 && true && frame_it==1)
        it = it + 1;

        oldK = IMG.K;
        IMG.SP_changed(:) = true;
        [IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy, newE] = split_move(IMG,1);
        E(end+1) = newE;
        converged = IMG.K - oldK < 2;

        if (dispOn)
            sfigure(1);
            subplot(1,1,1);
            im = zeros(size(oim));
            imagesc(IMG.label);
            title([num2str(it) ' - ' num2str(numel(unique(IMG.label)))]);

            sfigure(2);
            subplot(1,1,1);
            im = zeros(size(oim,1)+2*IMG.w, size(oim,2)+2*IMG.w, 3);
            im(IMG.w+1:end-IMG.w, IMG.w+1:end-IMG.w, :) = double(oim)/255;
            borders = is_border_valsIMPORT(double(reshape(IMG.label+1, [IMG.xdim IMG.ydim])));
            im = setPixelColors(im, find(borders), [1 0 0]);
            image(im,'parent',gca);
            drawnow;
        end
    end

    it = 0;
    converged = false;
    if (frame_it>1)
        IMG.SP_changed(:) = true;
        [IMG.K, IMG.label, IMG.SP, ~, IMG.max_UID, ~, ~, ~] = merge_move(IMG,1);
        [IMG.K, IMG.label, IMG.SP, ~, IMG.max_UID, ~, ~, ~] = split_move(IMG,10);
        [IMG.K, IMG.label, IMG.SP, ~, IMG.max_UID, ~, ~, ~] = switch_move(IMG,1);
        [IMG.K, IMG.label, IMG.SP, ~, IMG.max_UID, ~, ~, ~] = localonly_move(IMG,1000);
    end
    IMG.SP_changed(:) = true;
    IMG.alive_dead_changed = true;

    while (~converged && it<20)
        it = it + 1;
        times = zeros(1,5);

        if (~params.reestimateFlow)
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed1, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = localonly_move(IMG,1500);times(2)=toc;
            SP_changed0 = SP_changed1;
        else
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed0, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = local_move(IMG,1000);times(1)=toc;
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed1, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = localonly_move(IMG,500);times(2)=toc;
        end
        if (frame_it>1 && it<5)
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed2, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = merge_move(IMG,1);times(3)=toc;
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed3, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = split_move(IMG,1);times(4)=toc;
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed4, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = switch_move(IMG,1);times(5)=toc;
            IMG.SP_changed = SP_changed0 | SP_changed1 | SP_changed2 | SP_changed3 | SP_changed4;
        else
            IMG.SP_changed = SP_changed0 | SP_changed1;
        end

        E(end+1) = newE;
        converged = ~any(~arrayfun(@(x)(isempty(x{1})), {IMG.SP(:).N}) & IMG.SP_changed(1:IMG.K));

        if (dispOn)
            sfigure(1);
            im = zeros(size(oim));
            imagesc(IMG.label);
            title([num2str(it) ' - ' num2str(numel(unique(IMG.label)))]);

            sfigure(2);
            im = zeros(size(oim,1)+2*IMG.w, size(oim,2)+2*IMG.w, 3);
            im(IMG.w+1:end-IMG.w, IMG.w+1:end-IMG.w, :) = double(oim)/255;
            borders = is_border_valsIMPORT(double(reshape(IMG.label+1, [IMG.xdim IMG.ydim])));
            im = setPixelColors(im, find(borders), [1 0 0]);
            image(im,'parent',gca);
            drawnow;
        end
    end

    SP_UID = {IMG.SP(:).UID};
    mask = arrayfun(@(x)(isempty(x{1})), SP_UID);
    for m = find(mask)
        SP_UID{m} = -1;
    end
    sp_labels(:,:,frame_it) = reshape([SP_UID{IMG.label(IMG.w+1:end-IMG.w,IMG.w+1:end-IMG.w) +1}], size(oim,1), size(oim,2));

end



