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
function [Df, k0] = func_Df_k0 (g, Nodes, aggregateConfiguration, monomerRadius, boxLow, boxHigh, varargin)
    p = inputParser;
    %{
    %}
    addRequired(p, 'aggregateConfiguration', @(x) (size(x,2) == 3 || size(x,2) == 3));
    addRequired(p, 'monomerRadius', @(x) (x > 0) && isnumeric(x) && isscalar(x));
    addRequired(p, 'boxLow', @(x) (x > 0) && isnumeric(x) && isscalar(x));
    addRequired(p, 'boxHigh', @(x) (x > 0) && isnumeric(x) && isscalar(x));
    
    addParameter(p,'iterations',500);
    addParameter(p,'maxOutOfRadius',3); 
    addParameter(p,'outputcsv','');
    addParameter(p,'rngSeed',0);
    addParameter(p,'boxDims',0);
   
    parse(p, aggregateConfiguration, monomerRadius, boxLow, boxHigh, varargin{:});
    
   
    %% Generate data
    rng(p.Results.rngSeed,'twister'); % seed random number generator
    randNode = randi([1 length(aggregateConfiguration)],1,p.Results.iterations); % generated random integers

    % data arrays
    N_monomers = [];
    r_g = [];

    % box 
    %boxDims = [boxSize boxSize boxSize] * monomerRadius;

    for i = 1:length(randNode) % iterate over all random starting nodes
        % NOTE: replace print functions with waitbar()?
        %fprintf("Iteration: %d\n  Starting node: %d\n", i, randNode(i));
        %fprintf("  Search box: (%d, %d, %d)\n", boxDims(1), boxDims(2), boxDims(3));
        clear subAgg
        
        boxDims = randomBox(boxLow,boxHigh,monomerRadius);
        
        b = bfsearch(g, randNode(i));

        subAgg = Nodes(:, randNode(i)); % [x,y,z] of sub-aggregate compatible with Divjyot's RoG function
        outOfRadiusCounter = 0;
        
        for j = 1:size(b) % iterate over searched nodes
           if outOfRadiusCounter > p.Results.maxOutOfRadius
               break % note: this termination condition relies on order of BFS
           end

           if nodeIsInBox(Nodes, boxDims, randNode(i), b(j)) 
                subAgg = [subAgg Nodes(:, b(j))]; % probably should be replaced with a list
           else
               outOfRadiusCounter = outOfRadiusCounter + 1;
           end
       end

        subAgg = subAgg'; % Transpose why? i don't know
        new_n = size(subAgg,1)-1; 
        new_r = RoG(subAgg, monomerRadius);

        N_monomers = [N_monomers new_n];
        r_g = [r_g new_r];

        %fprintf("  N_monomers: %d\n  Radius of gyration: %.2f\n", new_n, new_r);
    end
    logrg = log(r_g) - log(monomerRadius);
    logN = log(N_monomers);
    coefficients = polyfit(logrg, logN, 1);
    
    Df = coefficients(1); 
    k0 = exp(coefficients(2));
end

%% Helper functions
function d = distance(a, b)
    temp = a-b;
    d = sqrt(dot(temp, temp));
end

function out = nodeIsInBox(Nodes, boxDims, center, candidate)
% center, candidate are indices
% boxDims is array containing x, y, z lengths of box
    temp = abs((Nodes(:,candidate)) - (Nodes(:,center)));
    for i = 1:3
        if temp(i) > boxDims(i)
            out = false;
            return
        end
    end
    out = true;
end

function out = randomBox(low, high, radius)
    out = randi([low high],1,3) * radius;
end