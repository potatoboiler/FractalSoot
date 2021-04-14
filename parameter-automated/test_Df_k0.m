%% Initialisation
delete(gcp('nocreate'));
clearvars; 
clf('reset');
% clc clears command line
margin = 2; %[nm]

%% "Hyperparameters"
% MAKE SURE FILENAME IS OF THE FORM "npart-nsamp-k0-df-iccmod.txt"
inputFile= ""; % Path to the input file
outputcsv = ''; % FILEPATH FOR CSV FILE
filetype = 'makowsk';
pool = parpool("threads");

%% Notes
% 120-180-1.3-1.8-1 --> optimal bounds are roughly [4, 7]
%

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
aggCon = monomersData(:,2:4); % aggregate configuration
%% Graph generation
% can also cache distance here as adjacency matrix
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
%% Conducting tests 
maxLength = floor(maxDistance(aggCon)/1);
minLength = 2;
% create test-space
test_space = 10*ones(maxLength, maxLength);

for high = minLength : maxLength
    tic
    parfor low = minLength : high
        fprintf("low: %d, high: %d\n", low, high);
        [Df, k0] = func_Df_k0(Nodes, BFStable, radius, low, high);
        x = MSE([Df k0], [1.8 1.3]);
        test_space(low, high) = min(10, x);
    end
    toc
end

mesh(test_space)
maxDistance(aggCon)
%% Helper functions
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