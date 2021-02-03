function [S, Obs] = generateSafeRegion_IRIS(obstacle_centers, safe_area_size, Gd_sensor, Fd_sensor)
%% Use Iris to find Safe regions
% Unit Cube:
verts = [   1.0 1.0 1.0;...
            -1.0 1.0 1.0;...
            1.0 -1.0 1.0;...
            -1.0 -1.0 1.0;...
            1.0 1.0 -1.0;...
            -1.0 1.0 -1.0;...
            1.0 -1.0 -1.0;...
            -1.0 -1.0 -1.0];
% Create obstacles:
obstacle_verts = zeros(3,8,size(obstacle_centers,1));
Obs = Polyhedron;
for i = 1:size(obstacle_centers,1)
    obstacle_verts(:,:,i) = (safe_area_size*verts + obstacle_centers(i,:))';
    Obs(i) = Polyhedron('V', obstacle_verts(:,:,i)');
end
% Choose the initial IRIS seed points
seeds = [[1.0; 0.0; 0.0], [0.0; 1.0; 0.0], [0.0; 0.0; 2.0], [0.0; 0.0; -2.0], [-2.0; -2.0; 0.0], [2.0; 2.0; 0.0], [0.75; 0.5; 0.0]];
% Build the initial IRIS regions
safe_regions = struct('ObsA', {}, 'Obsb', {}, 'point', {}); %, 'C', {}, 'd', {}, 'point', {});
for j = 1:size(seeds,2)
    s = seeds(:,j);
	[ObsA, Obsb,~,~,~] = iris.inflate_region(obstacle_verts, Gd_sensor(:,1:3), Fd_sensor, s, struct('require_containment', true, 'error_on_infeas_start', true));
	safe_regions(end+1) = struct('ObsA', ObsA, 'Obsb', Obsb, 'point', s); % 'C', C, 'd', d, 'point', s);
end

for i = 1:length(safe_regions)
    S(i) = Polyhedron('H',[safe_regions(i).ObsA safe_regions(i).Obsb]);
end