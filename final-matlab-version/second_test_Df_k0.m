% get list of filepaths
% subfigures each thingy, display the parameters for each
clear all;
warning('off','all');
clf('reset');
% clc clears command line
margin = 2; %[nm]

%% "Hyperparameters"
% MAKE SURE FILENAME IS OF THE FORM "npart-nsamp-k0-df-iccmod.txt"
inputFile= "C:\Users\Laurence\Documents\MATLAB\new-tests\120-180-1.3-1.8-0.txt"; % Path to the input file
outputcsv = ''; % FILEPATH FOR CSV FILE
filetype = 'fracmap';
low = 5;

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
toc

g = graph(edgeTail, edgeHead);

g.Nodes.X = aggCon(:,1);
g.Nodes.Y = aggCon(:,2);
g.Nodes.Z = aggCon(:,3);

Nodes = table2array(g.Nodes);
tic
BFStable = zeros(length(Nodes), length(Nodes));
for i = 1 : length(Nodes)
    BFStable(i,:) = bfsearch(g, i);
end
toc

rad = RoG(Nodes);

high = ceil(maxDistance(Nodes));
[Df, k0, errors] = func_Df_k0(Nodes, BFStable, radius, 7, ceil(rad));
maxDistance(Nodes);


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