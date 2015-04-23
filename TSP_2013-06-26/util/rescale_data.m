function [data] = rescale_data(data, p_Sigma, a_Sigma)

data(1,:) = data(1,:) / sqrt(p_Sigma(1));
data(2,:) = data(2,:) / sqrt(p_Sigma(2));
data(3,:) = data(3,:) / sqrt(a_Sigma(1));
data(4,:) = data(4,:) / sqrt(a_Sigma(2));
data(5,:) = data(5,:) / sqrt(a_Sigma(3));