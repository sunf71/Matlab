% The following variable links to the GSL libraries.
% On linux, this is probably not needed.
if (ismac)
    GSL_DIR = ' -I/usr/local/include/ -L/usr/local/lib/ ';
else
    GSL_DIR = '';
end

% TSP MEX functions

cd mex;

mex -O is_border_valsIMPORT.cpp
mex -O random_initIMPORT.cpp
mex -O populate_indices.cpp
mex -O SP_prop_init.cpp

eval(['mex ' GSL_DIR '-O localonly_move.cpp IMG.cpp NormalD.cpp SP.cpp -lgsl -lgslcblas -lm -lpthread']);
eval(['mex ' GSL_DIR '-O local_move.cpp IMG.cpp NormalD.cpp SP.cpp -lgsl -lgslcblas -lm -lpthread']);
eval(['mex ' GSL_DIR '-O switch_move.cpp IMG.cpp NormalD.cpp SP.cpp -lgsl -lgslcblas -lm -lpthread']);
eval(['mex ' GSL_DIR '-O merge_move.cpp IMG.cpp NormalD.cpp SP.cpp -lgsl -lgslcblas -lm -lpthread']);
eval(['mex ' GSL_DIR '-O split_move.cpp IMG.cpp NormalD.cpp SP.cpp -lgsl -lgslcblas -lm -lpthread']);

cd ..;

% optical flow

cd optical_flow_celiu/mex/
mex -O Coarse2FineTwoFrames.cpp OpticalFlow.cpp GaussianPyramid.cpp
cd ../../
