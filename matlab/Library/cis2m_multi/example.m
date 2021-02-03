clear
clc

A = [1 0 1; 0 0 1; 0 1 1];
B = [1 0; 1 1; 0 0];

n = 3; k = 10;

D = Polyhedron('H',[0 1 1]);
while (~D.isBounded && ~D.isEmptySet)
    G = -1 + (1+1)*rand(k,n);
    F = rand(k,1);
    D = Polyhedron('H',[G F]);
end

Gu = [eye(2); -eye(2)];
Fu = [ones(2,1)*0.01; ones(2,1)*0.03];
U = Polyhedron('H',[Gu Fu]);

L = 5;
cisMat = computeCIS(A,B,G,F,Gu,Fu,'ACC21b',L);
cis = Polyhedron('H', cisMat);


% Check if result is empty:
if (cis.isEmptySet)
    warning('result is empty');
end
% Check if result is numerically invariant:
if (~isInvariant(cis,A,B))
    warning('result not numerically invariant');
end
% Check if result is contained by the safe set:
if (~(cis <= D))
    warning('out of safe set');
end

volumeCIS = cis.volume;

%%%%%%%%%%%%%%%%%%%%%%%% Compute MCIS using MPT3 %%%%%%%%%%%%%%%%%%%%%%%%%
system = LTISystem('A',A,'B',B);
mcis = system.invariantSet('X',D,'U',U,'maxIterations',200);
volumeMCIS = mcis.volume;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(volumeMCIS)
disp(volumeCIS)

figure;
plot(D, 'color', 'blue', mcis, 'color', 'lightgray', cis, 'color', 'green')


cisMatImpl = computeImplicitCIS(A,B,G,F,Gu,Fu,'ACC21b',L);
cisImpl = Polyhedron('H', cisMatImpl);