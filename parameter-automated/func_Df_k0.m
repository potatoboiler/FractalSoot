%{
Input:
    aggregateConfiguration - 3xN array of N monomers
    monomerRadius - real number
    boxSize - real number, preferably integer
    varargin - 
        iterations - number of samples
        maxOutOfRadius - how many times a sample 
        outputcsv - not yet used
        rngSeed - 
%}
% TODO : investigate caching of distance info and pass to func
% NOTE : cache
function [Df, k0] = func_Df_k0 (g, Nodes, BFStable, monomerRadius, boxLow, boxHigh, varargin)
    p = inputParser;
    %{
        
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

    % data arrays
    % cannot pre-allocate, otherwise runtime increases from 0.2s to 5.0 seconds??
    N_monomers_array = [];
    radius_gyration_array = [];

    agg_size = length(Nodes);
    
    for i = 1:iterations % iterate over all random starting nodes
        % NOTE: replace print functions with waitbar()?
        % fprintf("Iteration: %d\n  Starting node: %d\n", i, randNode(i));
        % fprintf("  Search box: (%d, %d, %d)\n", boxDims(1), boxDims(2), boxDims(3));
        % clear subAgg % this can be safely removed
        
        %{  
            POSSIBILITY FOR INVESTIGATION:
            each node be its own object, containing list of : (other_node index, other_node
            distance)
            iterate over sorted list by distances
        %}
        boxDims = randomBox(boxLow,boxHigh,monomerRadius);
        
        b = BFStable(randNode(i),:);

        subAgg = zeros(agg_size, 3);
        subAgg(1,:) = Nodes(randNode(i),:); % [x,y,z] of sub-aggregate compatible with Divjyot's RoG function
        
        outOfRadiusCounter = 0;
        
        monomer_count = 1;
        
        for j = 1:size(b) % iterate over searched nodes
           if outOfRadiusCounter > maxOutOfRadius
               break % note: this termination condition relies on order of BFS
           end

           if nodeIsInBox(boxDims, Nodes(randNode(i),:), Nodes(b(j),:))
                monomer_count = monomer_count + 1;
                subAgg(monomer_count, :) = Nodes(b(j),:); % probably should be replaced with a list
           else
               outOfRadiusCounter = outOfRadiusCounter + 1;
           end
        end

        subAgg = subAgg(1:monomer_count,:);
        radius_gyration = RoG(subAgg, monomerRadius);

        N_monomers_array = [N_monomers_array monomer_count];
        radius_gyration_array = [radius_gyration_array radius_gyration];
        %fprintf("  N_monomers: %d\n  Radius of gyration: %.2f\n", new_n, new_r);
    end
    logrg = log(radius_gyration_array) - log(monomerRadius);
    logN = log(N_monomers_array);
    coefficients = polyfit(logrg, logN, 1);
    
    Df = coefficients(1); 
    k0 = exp(coefficients(2));
end

%% Helper functions
function d = distance(a, b)
    temp = a-b;
    d = sqrt(dot(temp, temp));
end

function out = nodeIsInBox(boxDims, center, candidate)
% center, candidate are subarrays of Nodes
% boxDims is array containing x, y, z lengths of box
% implements "manhattan distance" kinda
    temp = abs(candidate - center) - boxDims;
    out = ~any(temp(:) > 0);
end

function out = randomBox(low, high, radius)
    out = randi([low high],1,3) * radius;
end