function [] = compute_of(oim0, oim1, outname)

% add the optical flow path
addpath('optical_flow_celiu/');
addpath('optical_flow_celiu/mex/');

overwrite = false;

if (overwrite || ~exist(outname, 'file'))
    [~, flow.bvx, flow.bvy] = compute_flow(oim1, oim0);
    [~, flow.fvx, flow.fvy] = compute_flow(oim0, oim1);
    save(outname, 'flow');
end
