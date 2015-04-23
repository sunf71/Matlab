function [output] = setPixelColors(input, indices, color)
%function [output] = setPixelColors(input, indices, color)

N = size(input,1)*size(input,2);

output = input;
output(indices) = color(1);
output(indices + N) = color(2);
output(indices + 2*N) = color(3);