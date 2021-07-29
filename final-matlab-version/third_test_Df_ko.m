%% Initialisation
delete(gcp('nocreate'));
clearvars; 
clf('reset');
% clc clears command line
margin = 2; %[nm]

%{
keep min at 6
tests
%}  
%% "Hyperparameters"
% MAKE SURE FILENAME IS OF THE FORM "npart-nsamp-k0-df-iccmod.txt"
inputFile= "C:\Users\gamer\Documents\MATLAB\new-tests\new-fracs\200-300-0.68-2.29.txt"; % Path to the input file
outputcsv = ''; % FILEPATH FOR CSV FILE
filetype = 'fracmap';
% pool = parpool("threads");

%% Notes
% 120-180-1.3-1.8-1 --> optimal bounds are roughly [4, 7]
% 120 180 1.3 1.8 0 --> 
% 500 750 1.3 1.8 0 --> maxlen = 96; [ (30,50), (30,46) ]
% 1000 1500 1.3 1.8 0 --> maxlen = 133; [ (10,20), (15 --> 10) ]  [ (95,100), (95,100) ]
% 400 600 1.3 1.8 0 --> maxlen = 80; [30, 30] + small neighbourhood
% 200 300 1.3 1.8 0 --> maxlen = 48; [ (25,50), (25 + small neighbourhood) ]
% 150 225 1.3 1.8 0 --> 
% 900 1350 1.3 1.8 0 --> maxlen = 125. [neihgbourhood around (60,60), but not exactly there, roughly 2 radius]
% 
%% Extracting Data from file
tic
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
aggCon = monomersData(:,2:4); % aggregate configuration
toc
%% Graph generation
% can also cache distance here as adjacency matrix
tic
edgeTail = [];
edgeHead = [];
for i = 1:length(aggCon)
    for j = i+1:length(aggCon)
       if distance(aggCon(i,:), aggCon(j,:)) <= 2.1*radius
           edgeTail(end+1) = i;
           edgeHead(end+1) = j;
        end
    end
end

g = graph(edgeTail, edgeHead);

g.Nodes.X = aggCon(:,1);
g.Nodes.Y = aggCon(:,2);
g.Nodes.Z = aggCon(:,3);

Nodes = table2array(g.Nodes);

BFStable = zeros(length(Nodes), length(Nodes));
for i = 1 : length(Nodes)
    BFStable(i,:) = bfsearch(g, i);
end
toc
%% Conducting tests 
maxLength = floor(maxDistance(aggCon)/1);
minLength = 5;
% create test-space
test_space = 10*ones(maxLength, maxLength);
para_space = zeros(maxLength,maxLength,2);

stride = 100;
rad = ceil(RoG(Nodes));
%{
tic
for high = 2 : maxLength
    parfor low = minLength : high
        fprintf("low: %d, high: %d\n", low, high);
        [Df, k0] = func_Df_k0(Nodes, BFStable, radius, low, high);
        para_space(low, high, :) = [Df, k0];
        x = min(abs(10*(Df-1.84)), 10);
        test_space(low, high) = x; % test_space = errors
        %clf('reset');
    end
end
toc
%}

offset = 3;
high = rad;
x = offset*ones(high,1);
y = ones(high,1);
for low = 7 : 7
        fprintf("low: %d, high: %d\n", low, high);
        [Df, k0] = func_Df_k0(Nodes, BFStable, radius, low, high);
        para_space(low, high, :) = [Df, k0];
        x(low) = min([abs(Df-2.29) offset]);
        y(low) = Df;
        %clf('reset');
end
toc
hold on;
%plot(x);
%}


%mesh(test_space);
maxDistance(aggCon);
%% Helper functions
% TODO: normalize each value and _then_ sum+square
% TODO: take square root
function error = MSE(generated_params, sampled_params) 
    % both inputs should be pairs in the form [Dimension, prefactor] or
    % vice versa
    temp = generated_params - sampled_params;
    temp = 10*temp;
    error = dot(temp, temp);
end

function d = distance(a, b)
    temp = a-b;
    d = sqrt(dot(temp, temp));
end

function [d, node1, node2] = maxDistance(aggCon) 
    d = 0;
    for i = 1:length(aggCon) 
        for j = i : length(aggCon)
            a = distance(aggCon(i,:), aggCon(j,:)); 
            if a > d
                d = a;
                node1 = i;
                node2 = j;
            end
        end
    end
end