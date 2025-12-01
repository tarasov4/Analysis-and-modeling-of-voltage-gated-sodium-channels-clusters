%% set up workspace
clear
clc

working_path = '/Users/maddyammon/Documents/grad school work/sai lab rotation/minflux data/untreated/';
contents = dir([working_path, '*.tif']);

%import pixel size information
file_name = 'untreated.xlsx'; 
datainfo = readtable([working_path, file_name]); 
pixelsize = table2array(datainfo(:,3));
filename = datainfo(:,1);
Folder = pwd;

nbits = 8; % bit depth of input images
cutoff_multiplier = 7; % Clusters with area > median + cutoff_multiplier * std will be discarded

nfiles = length(contents); % # of files
cluster_mass_results = cell(nfiles,1); % Declare array to collect results
cluster_mass_values = [];
cluster_area_results = cell(nfiles,1); % Declare array to collect results
cluster_area_values = [];
cluster_density_results = cell(nfiles,1); % Declare array to collect results
cluster_density_values = [];
max_intensity = 2^nbits - 1;
xmax_mass = 50;% X axis limit for histograms
xmax_area = 10000;
xmax_density = .1;
bin_width = 1; % Bin width for histograms
bin_width_density = .001;
cluster_per_area = zeros(nfiles,1);

%% calculate cluster sizes
for counter1 = 1 : nfiles % Loop through files
    fn = contents(counter1).name; % Current file name
    cal = pixelsize(counter1); % Current calibration factor (nm)
    calibration = cal^2; % Area calibration factor (nm^2)
    im1 = imread([working_path, fn]); % Read image
    im1 = im1(:, :, 1); % Pull out only the red layer
    im1_binary = or(or(logical(im1 == 85), logical(im1 == 170)), logical(im1 == 255)); % Pixels with intensity values of 85, 170, 255 are considered signal-positive.
    im1_counts = zeros(size(im1));
    im1_counts(im1 == 85) = 1;
    im1_counts(im1 == 170) = 2;
    im1_counts(im1 == 255) = 3;

    cc = bwconncomp(im1_binary); % Perform KNN clustering

    nclusters = cc.NumObjects; % # of clusters
    pixels = cc.PixelIdxList; % pixel lists for clusters
    % Declare array to collect results
    cluster_mass = zeros(nclusters, 1); % Cluster Mass - total number of localizations in current clusters
    cluster_area = zeros(nclusters, 1); % Cluster Area - number of pixels in cluster * calibration factor
    % cluster_density = zeros(nclusters, 1); % Cluster Density - cluster mass normalized to cluster area
    for counter2 = 1 : nclusters % Loop through clusters
        pixel_values = im1_counts(pixels{counter2}); % Isolate, normalize values of pixels in current clusters
        cluster_mass(counter2) = sum(pixel_values); % Calculate mass of current cluster
        cluster_area(counter2) = numel(pixel_values) * calibration; % Cluster area
    end
    cluster_density = cluster_mass ./ cluster_area;

    cluster_mass_results{counter1} = cluster_mass; % store current image results in cell array
    cluster_mass_values = [cluster_mass_values; cluster_mass]; %#ok<AGROW>
    cluster_area_results{counter1} = cluster_area; % store current image results in cell array
    cluster_area_values = [cluster_area_values; cluster_area]; %#ok<AGROW>
    cluster_density_results{counter1} = cluster_density; % store current image results in cell array
    cluster_density_values = [cluster_density_values; cluster_density]; %#ok<AGROW>

end

% Remove statistical outliers
cutoff_mass = median(cluster_mass_values) + cutoff_multiplier * std(cluster_mass_values); % Max acceptable cluster mass
cutoff_area = median(cluster_area_values) + cutoff_multiplier * std(cluster_area_values); % Max acceptable cluster mass
cutoff_density = median(cluster_density_values) + cutoff_multiplier * std(cluster_density_values); % Max acceptable cluster mass
to_include_mass = logical(cluster_mass_values <= cutoff_mass);
to_include_area = logical(cluster_area_values <= cutoff_area);
to_include_density = logical(cluster_density_values <= cutoff_density);

cluster_mass_values = cluster_mass_values(to_include_mass); % Remove outliers from master list
cluster_area_values = cluster_area_values(to_include_area); % Remove outliers from master list
cluster_density_values = cluster_density_values(to_include_density); % Remove outliers from master list

% set up histogram bins
bin_edges_mass = 0 : bin_width : xmax_mass;
bin_edges_area = 0 : bin_width : xmax_area; 
bin_edges_density = 0 : bin_width_density : xmax_density;
[bin_counts_mass, bin_edges_mass] = histcounts(cluster_mass_values, bin_edges_mass, 'Normalization','cdf');
[bin_counts_area, bin_edges_area] = histcounts(cluster_area_values, bin_edges_area, 'Normalization','cdf');
[bin_counts_density, bin_edges_density] = histcounts(cluster_density_values, bin_edges_density, 'Normalization','cdf');
bin_centers_mass = bin_edges_mass(1 : end-1) + bin_width/2;
bin_centers_area = bin_edges_area(1 : end-1) + bin_width/2;
bin_centers_density = bin_edges_density(1 : end-1) + bin_width_density/2;

%% cdf plots
figure(1)
subplot(3,1,1)
plot(bin_centers_mass, bin_counts_mass)
set(gca, 'xlim', [0 xmax_mass])
xlabel('Cluster Mass [# localizations]')
ylabel('Cumulative Probability [0-1]')
title('Cluster Mass');

subplot(3,1,2)
plot(bin_centers_area, bin_counts_area)
set(gca, 'xlim', [0 xmax_area])
xlabel('Cluster Area [nm^2]')
ylabel('Cumulative Probability [0-1]')
title('Cluster Area');

subplot(3,1,3)
plot(bin_centers_density, bin_counts_density)
set(gca, 'xlim', [0 xmax_density])
xlabel('Cluster Density [# localizations / nm^2]')
ylabel('Cumulative Probability [0-1]')
title('Cluster Density');
exportgraphics(gcf, 'cluster_mass_CDF.pdf', 'ContentType','vector');

%% running k-means for mass (k=2)

idxk2 = kmeans(cluster_mass_values, 2); %run k-means and get a vector of values
indexk2_1 = logical(idxk2==1); %find all of the group 1 indices
indexk2_2 = logical(idxk2==2); %find all of the group 2 indices

[bin_counts_mass_all, ~] = histcounts(cluster_mass_values, bin_edges_mass, 'Normalization','cdf');

k1_max_val = max(cluster_mass_values(indexk2_1)); % Max cluster mass in k_1 group
k2_max_val = max(cluster_mass_values(indexk2_2)); % Max cluster mass in k_2 group

x_axis_breaks = sort([k1_max_val, k2_max_val]);
x_axis_break_k1 = find(bin_centers_mass <= x_axis_breaks(1), 1, 'last');
% x_axis_break_k2 = find(bin_centers_mass <= x_axis_breaks(2), 1, 'last');

%plot histogram for k means
figure(2)
histogram(cluster_mass_values(indexk2_1), 'BinEdges', bin_edges_mass);
hold on;
histogram(cluster_mass_values(indexk2_2), 'BinEdges', bin_edges_mass);
title('Cluster Mass, Kmeans(k = 2)')
xlabel('Cluster Mass [# localizations]')
hold off;
exportgraphics(gcf, 'cluster_mass_histogram.pdf', 'ContentType','vector');

figure(3)
plot(bin_centers_mass(1 : x_axis_break_k1), bin_counts_mass_all(1 : x_axis_break_k1), '-ok')
hold on
plot(bin_centers_mass(x_axis_break_k1 : end), bin_counts_mass_all(x_axis_break_k1 : end), '-or')
hold off
title('Cluster Mass, Kmeans(k = 2)')
xlabel('Cluster Mass [# localizations]')
ylabel('Cumulative Probability [0-1]')
exportgraphics(gcf, 'cluster_mass_CDF_k_means_k2.pdf', 'ContentType','vector');


%% running k-means for mass (k=3)

idxk3 = kmeans(cluster_mass_values, 3);

indexk3_1 = logical(idxk3==1); %find all of the group 1 indices
indexk3_2 = logical(idxk3==2); %find all of the group 2 indices
indexk3_3 = logical(idxk3==3); %find all of the group 3 indices

[bin_counts_mass_all, ~] = histcounts(cluster_mass_values, bin_edges_mass, 'Normalization','cdf');

k1_max_val = max(cluster_mass_values(indexk3_1)); % Max cluster mass in k_1 group
k2_max_val = max(cluster_mass_values(indexk3_2)); % Max cluster mass in k_2 group
k3_max_val = max(cluster_mass_values(indexk3_3)); % Max cluster mass in k_3 group

x_axis_breaks = sort([k1_max_val, k2_max_val, k3_max_val]);
x_axis_break_k1 = find(bin_centers_mass <= x_axis_breaks(1), 1, 'last');
x_axis_break_k2 = find(bin_centers_mass <= x_axis_breaks(2), 1, 'last');

%plot histogram for k means
figure(4)
histogram(cluster_mass_values(indexk3_1), 'BinEdges', bin_edges_mass);
hold on;
histogram(cluster_mass_values(indexk3_2), 'BinEdges', bin_edges_mass);
histogram(cluster_mass_values(indexk3_3), 'BinEdges', bin_edges_mass);
title('Cluster Mass, Kmeans(k = 3)')
xlabel('Cluster Mass [# localizations]')
hold off;
exportgraphics(gcf, 'cluster_mass_histogram_kmeans_k3.pdf', 'ContentType','vector');


figure(5)
plot(bin_centers_mass(1 : x_axis_break_k1), bin_counts_mass_all(1 : x_axis_break_k1), '-ok')
hold on
plot(bin_centers_mass(x_axis_break_k1 : x_axis_break_k2), bin_counts_mass_all(x_axis_break_k1 : x_axis_break_k2), '-or')
plot(bin_centers_mass(x_axis_break_k2 : end), bin_counts_mass_all(x_axis_break_k2 : end), '-om')
hold off
title('Cluster Mass, Kmeans(k = 3)')
xlabel('Cluster Mass [# localizations]')
ylabel('Cumulative Probability [0-1]')
exportgraphics(gcf, 'cluster_mass_CDF_kmeans_k3.pdf', 'ContentType','vector');



%% running k-means for density (k=2)

idxk2_density = kmeans(cluster_density_values, 2); %run k-means and get a vector of values
indexk2_density_1 = logical(idxk2_density==1); %find all of the group 1 indices
indexk2_density_2 = logical(idxk2_density==2); %find all of the group 2 indices

[bin_counts_density_all, ~] = histcounts(cluster_density_values, bin_edges_density, 'Normalization','cdf');

k1_density_max_val = max(cluster_density_values(indexk2_density_1)); % Max cluster density in k_1 group
k2_density_max_val = max(cluster_density_values(indexk2_density_2)); % Max cluster density in k_2 group

x_axis_breaks_density = sort([k1_density_max_val, k2_density_max_val]);
x_axis_break_k1_density = find(bin_centers_density <= x_axis_breaks_density(1), 1, 'last');

%plot histogram for k means
figure(6)
histogram(cluster_density_values(indexk2_density_1), 'BinEdges', bin_edges_mass);
hold on;
histogram(cluster_density_values(indexk2_density_2), 'BinEdges', bin_edges_mass);
title('Cluster Density, Kmeans(k = 2)')
xlabel('Cluster Density [# localizations / nm^2]')
hold off;
exportgraphics(gcf, 'cluster_density_histogram_kmeans_k2.pdf', 'ContentType','vector');

figure(7)
plot(bin_centers_density(1 : x_axis_break_k1_density), bin_counts_density_all(1 : x_axis_break_k1_density), '-ok')
hold on
plot(bin_centers_density(x_axis_break_k1_density : end), bin_counts_density_all(x_axis_break_k1_density : end), '-or')
hold off
title('Cluster Density, Kmeans(k = 2)')
xlabel('Cluster Density [# localizations / nm^2]')
ylabel('Cumulative Probability [0-1]')
exportgraphics(gcf, 'cluster_density_CDF_kmeans_k2.pdf', 'ContentType','vector');


%% running k-means for density (k=3)

idxk3_density = kmeans(cluster_density_values, 3);

indexk3_density_1 = logical(idxk3_density==1); %find all of the group 1 indices
indexk3_density_2 = logical(idxk3_density==2); %find all of the group 2 indices
indexk3_density_3 = logical(idxk3_density==3); %find all of the group 3 indices

[bin_counts_density_all, ~] = histcounts(cluster_density_values, bin_edges_density, 'Normalization','cdf');

k1_density_max_val = max(cluster_density_values(indexk3_density_1)); % Max cluster mass in k_1 group
k2_density_max_val = max(cluster_density_values(indexk3_density_2)); % Max cluster mass in k_2 group
k3_density_max_val = max(cluster_density_values(indexk3_density_3)); % Max cluster mass in k_3 group

x_axis_breaks_density = sort([k1_density_max_val, k2_density_max_val, k3_density_max_val]);
x_axis_break_k1_density = find(bin_centers_density <= x_axis_breaks_density(1), 1, 'last');
x_axis_break_k2_density = find(bin_centers_density <= x_axis_breaks_density(2), 1, 'last');

%plot histogram for k means
figure(8)
histogram(cluster_density_values(indexk3_density_1), 'BinEdges', bin_edges_density);
hold on;
histogram(cluster_density_values(indexk3_density_2), 'BinEdges', bin_edges_density);
histogram(cluster_density_values(indexk3_density_3), 'BinEdges', bin_edges_density);
title('Cluster Density, Kmeans(k = 3)')
xlabel('Cluster Density [# localizations / nm^2]')
hold off;
exportgraphics(gcf, 'cluster_density_histogram_kmeans_k3.pdf', 'ContentType','vector');

figure(9)
plot(bin_centers_density(1 : x_axis_break_k1_density), bin_counts_density_all(1 : x_axis_break_k1_density), '-ok')
hold on
plot(bin_centers_density(x_axis_break_k1_density : x_axis_break_k2_density), bin_counts_density_all(x_axis_break_k1_density : x_axis_break_k2_density), '-or')
plot(bin_centers_density(x_axis_break_k2_density : end), bin_counts_density_all(x_axis_break_k2_density : end), '-om')
hold off
title('Cluster Density, Kmeans(k = 3)')
xlabel('Cluster Density [# localizations / nm^2]')
ylabel('Cumulative Probability [0-1]')
exportgraphics(gcf, 'cluster_density_CDF_kmeans_k3.pdf', 'ContentType','vector');

%% save results
save('results.mat', "cluster_mass_results", "cluster_mass_values", "cutoff_mass",...
    "cluster_area_results", "cluster_area_values", "cluster_density_results", "cutoff_area",...
    "cluster_density_values", "cutoff_density");

save('clusterPerArea.mat', "cluster_per_area");