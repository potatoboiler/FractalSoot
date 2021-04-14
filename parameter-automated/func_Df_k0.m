%{
Box-counting algorithm implemented by Laurence Lu - nickname "FracBox"
Last modified: 4/13/2021

This algorithm iterates over randomly selected monomers from the aggregate.
Taking a box around this central monomer, connectivity to the central monomer is determined by iterating through the breadth-first search of monomers radiating outwards from the central monomer.
Collecting the radius of gyration of the sub-aggregate inside the box and the number of monomers together, these values are plotted against each other.
The slope from the resulting plot is the fractal dimension, the intercept is the log of the prefactor.

NOTES:
This function requires a precomputed table of bread-first searches along an aggregate
Connectivity information must be provided by a wrapper function or otherwise.
The breadth-first search means that nodes closer to the central monomer in
the aggregate connectivity-wise will be traversed first, guaranteeing
faster runtime. The order guarantees that very few, if any at all, monomers
are excluded from the sub-aggregate.
%}
% TODO : investigate caching of distance info and pass to func
function [fractal_dimension, prefactor] = func_Df_k0 (Nodes, BFStable, monomerRadius, boxLow, boxHigh, varargin)
    p = inputParser;
    %{
        Input:
            Nodes - N*3 array of monomer coordinates
            BFStable - N*N array, indexed along axis 1, each element is a breadth-first search list
            monomerRadius - real number
            boxLow - lower bound on possible box sizes; real number, preferably integer
            bowHigh - upper bound on possible box sizes; real number, preferably integer
            varargin - 
                iterations - number of samples
                maxOutOfRadius - how many times the algorithm will check 
                outputcsv - not yet used
                rngSeed - seed for random number generator
    %}
    addRequired(p, 'Nodes', @(x) (size(x,2) == 3 || size(x,2) == 3));
    addRequired(p, 'monomerRadius', @(x) (x > 0) && isnumeric(x) && isscalar(x));
    addRequired(p, 'boxLow', @(x) (x > 0) && isnumeric(x) && isscalar(x));
    addRequired(p, 'boxHigh', @(x) (x > 0) && isnumeric(x) && isscalar(x));
    
    addParameter(p,'iterations',500);
    addParameter(p,'maxOutOfRadius',3); 
    addParameter(p,'outputcsv','');
    addParameter(p,'rngSeed',0);
    addParameter(p,'boxDims',0);
   
    parse(p, Nodes, monomerRadius, boxLow, boxHigh, varargin{:});
    
    maxOutOfRadius = p.Results.maxOutOfRadius;
    iterations = p.Results.iterations;
    
    %% Generate data
    rng(p.Results.rngSeed,'twister'); % seed random number generator
    randNode = randi([1 length(Nodes)],1,iterations); % generated random integers

    % cannot pre-allocate, otherwise runtime increases from 0.2s to 5.0 seconds??
    N_monomers_array = [];
    radius_gyration_array = [];

    agg_size = length(Nodes);
    
    for i = 1:iterations % iterate over all random starting nodes
        % NOTE: replace print functions with waitbar()?
        % fprintf("Iteration: %d\n  Starting node: %d\n", i, randNode(i));
        % fprintf("  Search box: (%d, %d, %d)\n", boxDims(1), boxDims(2), boxDims(3));
        % clear subAgg % this can be safely removed
                
        % generates 3 random integers (box dimensions), then scales by radius 
        % to guarantee that at least a few monomers will be contained in the box
        boxDims = randi([boxLow boxHigh],1,3) * radius;
        
        b = BFStable(randNode(i),:); % refer to notes about input arguments

        subAgg = zeros(agg_size, 3); % pre-allocate space for arrays
        subAgg(1,:) = Nodes(randNode(i),:); % assign [x y z] of sub-aggregate; 
        
        outOfRadiusCounter = 0; % sentinel value for traversing BFS
        
        monomer_count = 1; % number of monomers in sub-aggregate
        
        for j = 1:length(b) % iterate over BFS list
           if outOfRadiusCounter > maxOutOfRadius % early termination condition
               break
           end
           
            % abs(Nodes(randNode(i),:) - Nodes(b(j),:)) = absolute
            % difference in x,y,z components between monomers
            % If *any* components are greater than zero after subtracting
            % the dimension of the box, the second monomer lies outside of the box. 
            % This expression is negated (in hopes that shortcircuiting
            % via any() is slightly faster than all())
           if ~any( (abs(Nodes(randNode(i),:) - Nodes(b(j),:)) - boxDims) > 0)
                monomer_count = monomer_count + 1;
                subAgg(monomer_count,:) = Nodes(b(j),:);
           else
               outOfRadiusCounter = outOfRadiusCounter + 1;
           end
        end

        subAgg = subAgg(1:monomer_count,:); % truncate rows of zeros
        radius_gyration = RoG(subAgg, monomerRadius); % compute radius of gyration

        N_monomers_array = [N_monomers_array monomer_count]; % append monomer_count
        radius_gyration_array = [radius_gyration_array radius_gyration]; % append radius_gyration
        %fprintf("  N_monomers: %d\n  Radius of gyration: %.2f\n", new_n, new_r);
    end
    logrg = log(radius_gyration_array) - log(monomerRadius); % scale radius of gyration values to fit fractal law
    logN = log(N_monomers_array);
    coefficients = polyfit(logrg, logN, 1);
    
    fractal_dimension = coefficients(1); % Df
    prefactor = exp(coefficients(2)); % k0
end