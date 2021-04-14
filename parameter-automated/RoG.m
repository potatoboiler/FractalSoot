% Written by Divjyot Singh (modified by Laurence Lu)
% Last Modified: April 10, 2021
% This function calculates the radius of gyration for a soot aggregate

% 'aggregate' is an N*3 array which stores the x,y,z coordinates of
% the N monomers in a soot aggregate. 'radius' is the radius of the monomer
% of the monondisperse soot aggregate. All units are nm. All variables are
% float.

% the function returns 'radiusOfGyration', which is the radius of gyration
% of a given soot aggregate. For aggregates with N>1, this function assumes
% the monomers to be point mass. This approximation is not expected to
% significantly impact the results.

% modification: marginal speed improvement (roughly 10-30%)

function radiusOfGyration = RoG(aggregate, radius)
    Ns = size(aggregate,1);
    if Ns == 1
        radiusOfGyration = sqrt(3/5)*radius;
        return
    end

    center = mean(aggregate,1);

    x = aggregate(:, 1) - center(1);
    y = aggregate(:, 2) - center(2);
    z = aggregate(:, 3) - center(3);    
    Rt = dot(x,x) + dot(y,y) + dot(z,z); % equivalent to squaring each axis and summing
    radiusOfGyration = (Rt/Ns)^0.5; 
end

