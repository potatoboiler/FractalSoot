% Notes:
%  I strongly suspect this can be made more efficient
% TODO:
%  Figure out why the parameters dont seem to match?
%  Automate this? 
%% Initialisation

clearvars % clears existing variables in the workspace
clf('reset') % reset figure window

% MAKE SURE FILENAME IS OF THE FORM "npart-nsamp-k0-df-iccmod.txt"
inputFile= "C:\Users\gamer\Documents\MATLAB\convexity\fractal_properties\Makowski\compile\1000-1000-1.3-1.8-0.txt"; %"C:\Users\gamer\Documents\MATLAB\divjyotscode\FracMap_600.txt"; % Path to the input file
margin = 2; %[nm]

%% "Hyperparameters"
distanceCheck = 'box'; % 'box' or 'radius'
iterations = 100; % how many sub-aggregates are sampled
maxOutOfRadius = 3; % how many times distance is checked before giving up
filetype = 'makowsk'; % either 'makowsk' or 'fracmap'
rngSeed = 0;

box_low = 3;
box_high = 100;

rad_low = 5;
rad_high = 15;

%% Extracting Data from file

monomersFile = fopen (inputFile,'r');
if (monomersFile == -1)
    fprintf('Check the input file path\n');
    return;
end

monomersData = 0;
radius = 0;
aggCon = 0;

if(strcmp('fracmap', filetype))    
    monomersData = textscan(monomersFile,'%f %f %f %f %f %f');
    refractiveIndex = monomersData(:,[5 6]);
elseif(strcmp('makowsk',filetype))
    monomersData = textscan(monomersFile,'%f %f %f %f');
end
monomersData = cell2mat(monomersData);
radius = monomersData(1,1);
aggCon = monomersData(:,[2 3 4]); % aggregate configuration

%% Generate graph
edgeTail = [];
edgeHead = [];
for i = 1:length(aggCon)
    for j = i+1:length(aggCon)
       if distance(aggCon(i,:), aggCon(j,:)) <= 2.1*radius
           edgeTail(end+1) = i;
           edgeHead(end+1) = j;
           % [i j];
        end
    end
end

g = graph(edgeTail, edgeHead);
g.Nodes.X = aggCon(:,1);
g.Nodes.Y = aggCon(:,2);
g.Nodes.Z = aggCon(:,3);
%plot(g)

%% Generate data

% data arrays
N_monomers = [];
r_g = [];

% radial distance check 
rng(rngSeed,'twister'); % seed random number generator
randNode = randi([1 length(aggCon)],1,iterations); % generated random integers
randRadius = randi([rad_low rad_high],1,iterations) * radius; % used for 'radius' proximity checks

% box 
boxDims = randomBox(box_low,box_high,radius);

for i = 1:length(randNode) % iterate over all random starting nodes
    %fprintf("Iteration: %d\n  Starting node: %d\n  Search radius: %d\n", i, randNode(i), randRadius(i));
    fprintf("Iteration: %d\n  Starting node: %d\n", i, randNode(i));
    clear subAgg
    %  Search width/radius: %d\n
    
    b = bfsearch(g, randNode(i));
    subAgg = table2array(g.Nodes(randNode(i),:))'; % [x,y,z] of sub-aggregate compatible with Divjyot's RoG function
    outOfRadiusCounter = 0;
    
    if strcmp('box', distanceCheck)
        fprintf("  Search box: (%d, %d, %d)\n", boxDims(1), boxDims(2), boxDims(3));
    else
        fprintf("  Search box: %d\n", randRadius(i));
    end
    
    for j = 1:size(b) % iterate over searched nodes
       if outOfRadiusCounter > maxOutOfRadius
           break % note: this termination condition relies on order of BFS
       end
       
       if strcmp('box', distanceCheck)
           if nodeIsInBox(g, boxDims, randNode(i), b(j)) 
                subAgg = [subAgg table2array(g.Nodes(b(j),:))']; % probably should be replaced with a list
           else
               outOfRadiusCounter = outOfRadiusCounter + 1;
           end
       elseif strcmp('radius', distanceCheck)
           if nodeDistance(g, randNode(i), b(j)) < randRadius(i) % nodeDistance(g, randNode(i), b(j)) <= randRadius(i)
                subAgg = [subAgg table2array(g.Nodes(b(j),:))']; % probably should be replaced with a list
           else
               outOfRadiusCounter = outOfRadiusCounter + 1;
           end
       end
    end
    
    subAgg = subAgg'; % Transpose why? i don't know
    new_n = size(subAgg,1)-1; 
    new_r = RoG(subAgg, radius);
    
    N_monomers = [N_monomers new_n];
    r_g = [r_g new_r];
    
    fprintf("  N_monomers: %d\n  Radius of gyration: %.2f\n", new_n, new_r);
end

%% Data processing

logrg = log(r_g) - log(radius);
logN = log(N_monomers);

coefficients = polyfit(logrg, logN, 1);
xFit = linspace(min(logrg), max(logrg), 1000);
yFit = polyval(coefficients, xFit);

eq_str = sprintf('y = %3.3f * x + %3.3f', coefficients(1), coefficients(2))
fprintf("Fractal dimension: %f\nPrefactor: %f\n", (coefficients(1)), exp(coefficients(2)));

figure(5)
hold on
title('Log-log graph of number of monomers in a sub-aggregate against its radius of gyration')
xlabel('log( R_g/a )');
ylabel('log( N )');
scatter(logrg, logN)
plot(xFit, yFit)
text(2,4, eq_str, 'FontSize', 16)
hold off


%% Output data to file

[~, name, ~] = fileparts(inputFile);
fileparams = textscan(name, "%f-%f-%f-%f-%f");
params = [distanceCheck iterations maxOutOfRadius filetype rngSeed " " box_low box_high " " rad_low rad_high];
data = [exp(coefficients(1)) exp(coefficients(2))];
outputRow = [name " " fileparams " " data " " params];

writematrix(outputRow,'C:\Users\gamer\Documents\MATLAB\convexity\fractal_properties\data.csv','WriteMode','append');

fclose('all');
clear inputFile
clear monomersFile
clear monomersData

%% Functions
function d = distance(a, b)
    d = sqrt(dot(a-b, a-b));
end

function d = nodeDistance(g, a,b) 
% assume that a,b are node indices
% assume that g has X,Y,Z coordinates per each node
    temp = table2array(g.Nodes(a,:)) - table2array(g.Nodes(b,:));
    d = sqrt(dot(temp, temp));
end

function out = randomBox(low, high, radius)
    out = randi([low high],1,3) * radius;
end

function out = nodeIsInBox(g, boxDims, center, cand)
% center, cand are indices
% boxDims is array containing x, y, z lengths of box
    temp = abs(table2array(g.Nodes(cand,:)) - table2array(g.Nodes(center,:)));
    for i = 1:3
        if temp(i) > boxDims(1)
            out = false;
            return
        end
    end
    out = true;
end
