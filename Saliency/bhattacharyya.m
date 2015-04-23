% BHATTACHARYYA(histogram1, histogram2)
% compute the BHATTACHARYYA distance between 2 histograms
% where each histogram is a 1xN vector
% 
% Based on the estimation presented in 
% "Real-Time Tracking of Non-Rigid Objects using Mean Shift"
%  D. Comaniciu, V. Ramesh & P. Meer (IEEE CVPR 2000, 142-151)
%
% N.B. both histograms must be normalised
% (i.e. bin values lie in range 0->1 and SUM(bins(i)) = 1
%       for i = {histogram1, histogram2} )
%
% Author / Copyright: T. Breckon, October 2005.
% School of Informatics, University of Edinburgh.
% License: http://www.gnu.org/licenses/gpl.txt

function bdist = bhattacharyya(histogram1, histogram2)

    % get number of bins 
    % (i.e. dimension 2 from Nx1 inputs)

    bins = size(histogram1, 2);
    
    % estimate the bhattacharyya co-efficient
    
    bcoeff = 0;
    for i=1:bins

        bcoeff = bcoeff + sqrt(histogram1(i) * histogram2(i));

    end
    
    % get the distance between the two distributions as follows
    
    bdist = sqrt(1 - bcoeff);

end