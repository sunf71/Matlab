function [sdf] = calc_sdf(regions)
%CALC_SDF    Calculates signed distance function
%   CALC_SDF calculates a signed distance function. 'regions' is a binary
%   mask indicating the binary regions to infer an SDF from.  The SDF is
%   calculated for a narrow band indicated by band_size.  Additionally, if
%   we have an SDF already calculated for most of the region, and only part
%   of the region changes, we can accelerate the process a lot. In that
%   case, mask indicates all points that have changed, and oldsdf indicates
%   the previously calculated SDF.
%
%   [sdf] = calc_sdf(regions, band_size, mask, oldsdf)
%
%   Required Parameters:
%     regions - a matrix containing the labelmap of the region. Only the
%       signs of regions has to be correct (positive in one region, negative
%       in the other region)
%     band_size - the narrow band size used
%
%   Optional Parameters:
%     mask - binary mask indicating pixels that have changed or pixels that
%       are near ones that changes and need to be updated
%     oldsdf - the previously calculated SDF used with mask
%
%   Copyright 2011. MIT. All Rights Reserved.
%   Written by Jason Chang 09-06-2011

%   Dependent files:
%     calc_sdfIMPORT.cpp

regions = int32(regions);
sdf = calc_sdfIMPORT(regions);