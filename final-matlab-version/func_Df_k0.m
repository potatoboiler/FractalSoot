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
% TODO : save the error for the polynomial fits [plot each line regression from each scatter plot]
%       - different methods of variation: one plot with varying just top
%       range, or just bot range, or both
% TODO : randomly : average over different runs of the same plots
% TODO : self-implement BFS, recurse, change to radius?
function [fractal_dimension, prefactor, errors] = func_Df_k0 (Nodes, BFStable, monomerRadius, boxLow, boxHigh, fracDim, pref, varargin)
    p = inputParser;
    %{
        Input:
            Nodes - N*3 array of monomer coordinates
            BFStable - N*N array, indexed along axis 1, each element is a breadth-first search list
            monomerRadius - real number
            boxLow - lower bound on possible box sizes; real number, preferably integer
            bowHigh - upper bound on possible box sizes; real number, preferably integer
            fracDim - fractal dimension value to be compared against
            pref - prefactor value to be compared against
            varargin - 
                iterations - number of samples
                maxOutOfRadius - how many times the algorithm will check 
                outputcsv - not yet used
                rngSeed - seed for random number generator
    %}
    addRequired(p, 'Nodes', @(x) (size(x,2) == 3 || size(x,2) == 3));
    addRequired(p, 'BFStable');
    addRequired(p, 'monomerRadius', @(x) (x > 0) && isnumeric(x) && isscalar(x));
    addRequired(p, 'boxLow', @(x) (x > 0) && isnumeric(x) && isscalar(x));
    addRequired(p, 'boxHigh', @(x) (x > 0) && isnumeric(x) && isscalar(x));
    
    addParameter(p,'iterations',500);
    addParameter(p,'maxOutOfRadius',100); 
    addParameter(p,'outputcsv','');
    addParameter(p,'rngSeed',0);
    addParameter(p,'boxDims',0);
   
    parse(p, Nodes, BFStable, monomerRadius, boxLow, boxHigh, varargin{:});
    
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
        boxDims = randi([boxLow boxHigh],1,3); %* monomerRadius;
        %a = randi([boxLow boxHigh],1,3) * monomerRadius;
        %boxDims = [a a a];
        
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
           %temp = (Nodes(randNode(i),:) - Nodes(b(j),:));
           %if (temp'*temp) < boxDims(1)*boxDims(1)
                monomer_count = monomer_count + 1;
                subAgg(monomer_count,:) = Nodes(b(j),:);
           else
               outOfRadiusCounter = outOfRadiusCounter + 1; % consider that this might undersample
               % recursively search down, instead, because of this
               % undersampling between hitting one wall first
           end
        end

        subAgg = subAgg(1:monomer_count,:); % truncate rows of zeros
        radius_gyration = RoG(subAgg, monomerRadius); % compute radius of gyration

        N_monomers_array = [N_monomers_array monomer_count]; % append monomer_count
        radius_gyration_array = [radius_gyration_array radius_gyration]; % append radius_gyration
        %fprintf("  N_monomers: %d\n  Radius of gyration: %.2f\n", new_n, new_r);
    end

    logrg = log10(radius_gyration_array) - log(monomerRadius); % scale radius of gyration values to fit fractal law
    logN = log10(N_monomers_array);
    
    % polyfit() is faster and more efficient
    % we want errors, so we use regress() for convenience
    [coefficients, raw_errors] = regress(logN', cat(2, logrg', ones(length(logrg))));
    %coefficients = polyfit(logrg, logN, 1); raw_errors = [];

    
    fractal_dimension = coefficients(1); % Df
    %prefactor = 10.^(coefficients(2)); % k0
    prefactor = agg_size / (RoG(Nodes, monomerRadius)/monomerRadius).^fractal_dimension;
   
    
    
    raw_errors = (raw_errors(1:2, 2) - raw_errors(1:2,1))/2;
    errors = raw_errors;
    errors(2) = errors(2) * prefactor; % propagate error for exponentiated term
    
    %errors = [ Df_error, k0_error ]
    
    %% Plotting tools
    
    % These can be safely commented out. 
    % These lines plot the data gathered from the previous section.
    xFit = linspace(min(logrg), max(logrg), 1000);
    yFit = xFit*coefficients(1) + coefficients(2);
    
    fprintf("Fractal dimension: %f\nPrefactor: %f\n", fractal_dimension, prefactor);
    eq_str = sprintf('y = %3.3f * x + %3.3f', fractal_dimension, coefficients(2));
    
    %str_arr = ['./figures/all/', num2str(length(Nodes)),'-',num2str(iterations),'-', num2str(boxLow), '-', num2str(boxHigh),'.png'];
    str_arr = ['./figures/all/', num2str(length(Nodes)),'-', pref, '-', fracDim, '.png'];
    
    filename = join(str_arr);
    
    info_str = sprintf('Dimension: %.2f %c %.2f\nPrefactor: %.2f %c %.2f', fractal_dimension, char(177), errors(1),...
        prefactor, char(177), errors(2));
    
    
    figure(1);
    hold on;
    grid on;
    scatter(logrg, logN, '.', 'r')
    
    %scatter(logrg, logN,'r')
    
    idealLine = xFit*1.68 + log10(1.08);
    %plot(xFit, idealLine, 'b')
    title('Log-log graph of sub-aggregate monomer count vs. radius of gyration');
    xlabel('log( R_g/a )');
    ylabel('log( N )');
    plot(xFit, yFit,'m')
    %plot(xFit,yFit+raw_errors(2),'m--',xFit,yFit-raw_errors(2),'m--')
    %text(2,4, eq_str, 'FontSize', 16)
    annotation('textbox', [0.6, 0.1, 0.1, 0.1], 'String', eq_str)
    annotation('textbox', [0.6, 0.2, 0.1, 0.1], 'String', info_str)
    saveas(gcf,filename);
    hold off  
    %}
end
