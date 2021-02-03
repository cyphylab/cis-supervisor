function [Graph] = generateDSetsMap(dsets)
% Function to generate the CISs graph with extra information
%
% Outputs
% Graph: Structure with fields:
%   Index: Index of the Node
%   Region: Polytope of the represented region
%   Intersections: Cell array of the intersection with other regions
%   Adjacency: Vector of flags of the node adjacency information
%   Volumes: Vector with information regarding the intersection volume.
%
%
N = length(dsets);

AMatrix = zeros(N);
AdjacentVolumes = zeros(N);
Intersection = cell(N, N);

Graph = cell(N,1);
for i = 1:N
    Graph{i}.Index = i;
    Graph{i}.Region = dsets{i};
    
    for j = i + 1:length(dsets)
        Intersection{i, j} = dsets{i}.Polyhedron & dsets{j}.Polyhedron;
        
        acc_int = Intersection{i, j}.projection(7:9, 'ifourier');
        vol = acc_int.volume();
        if (~Intersection{i, j}.isEmptySet)
            AMatrix(i, j) = 1;
            AMatrix(j, i) = 1;
            AdjacentVolumes(i, j) = vol;
            AdjacentVolumes(j, i) = vol;
        end
    end
    Graph{i}.Intersections = Intersection(i, :);
    Graph{i}.Adjacency = AMatrix(i, :);
    Graph{i}.Volumes = AdjacentVolumes(i, :);
end


end