function weights = CalcNeighborWeight(radius,sigma)
        sum = 0;
        weights = zeros(2*radius+1,1);
        
        for i=-radius:radius
            dis = abs(i);
            weights(i+radius+1) = (1/sqrt(2*pi*sigma))*exp((-dis*dis)/(2*sigma*sigma));
            sum = weights(i+radius+1) + sum;
        end
        
        weights = weights./sum;
end