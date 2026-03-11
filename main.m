clc
clear

%% If you use our toolbox or if you find our research helpful, please cite the following paper (thanks for your support):
%% Chen HL, Shen YZ. and Zhang L. 2025. Iterative Weighted Least Squares Solution to InSAR Stratified Tropospheric Phase Delay. IEEE Journal of Selected Toptics in Applied Earth Observations and Remote Sensing, 18, 19594 - 19605
%% Any feedback is welcome! (2210929@tongji.edu.cn or yzshen@tongji.edu.cn)

%% 1. Prepare input data
addpath '.\weighted_InSAR strified_delay_estimation'
% ts_file='timeseries_SET.h5';% Sentinel time series derived from MintPy and ISCE2
% s=3;% Descending take as 3; Ascending take as 1
% [phase,longitude1,latitude1,height1]= read_H5data(ts_file,s);
% save InSAR_input.mat  phase longitude1 latitude1 height1
% 
% 

load InSAR_input.mat  phase longitude1 latitude1 height1

%% 2. Predefined parameters
 Number_segments=40;%Used for calculating covariance
 sample=10;% Sparse sampling, reducing the computational burden

%% 3. Set reference point
refpoint=[1200,520];% The reference point can be chosen arbitrarily. The results under different reference points vary very little which is also the advantage of this weighted method.
height1=height1-height1(refpoint(1,1),refpoint(1,2));
phase=phase-phase(refpoint(1,1),refpoint(1,2),:);
imagesc(height1);colorbar;colormap(jet)

%% 4. DEM normalization
height1=height1./(max(height1(:))-min(height1(:)));

%% 5. Set up the function model for estimating stratification
order=3; %  1:first-order linear;2:second-order linear;3:third-order linear(When there is a large height difference, adopt)

%% 6. Establish a differential time series as the observed values to weaken the influence of deformation.
for i=1:size(phase,3)-1
phase_diff(:,:,i)=phase(:,:,i+1)-phase(:,:,i);
end
clear phase

%% 7. Data downsampling
latitude1_sample = latitude1(1:sample:end, 1:sample:end);
longitude1_sample = longitude1(1:sample:end, 1:sample:end);
phase_diff_sample = phase_diff(1:sample:end, 1:sample:end,:);
height1_sample = height1(1:sample:end, 1:sample:end);
clear latitude1 longitude1 
location=[longitude1_sample(:) latitude1_sample(:)];
location(any(isnan(location),2),:) = [];
[X,Y]= Pub_BLH2XYZ(location(:,2),location(:,1));
 XY=[X Y];
 clear Y X Y1 X1
 DS=squareform(pdist(XY));
 DS=single(DS);
  
clear location longitude1_sample latitude1_sample
%% 8. Estimate the parameters of the linear model
for i=1:size(phase_diff_sample,3)
     l=phase_diff_sample(:,:,i);
     [par_weight(i,:),par_unweight(i,:),hgt]=Par_estimation(l,height1_sample,DS,Number_segments,order);
     i
 end

%% 9. time series correction
 Str_correction=calculate_correction_timeseries(par_weight,height1,order);
 imagesc( Str_correction(:,:,2));colorbar
%% 10. Outlooks
% For situations where stratification needs to be estimated, after the window has been divided, the model coefficients can be estimated using weighted methods.
% The remaining severe turbulence has not been dealt with yet. We will update the code for handling turbulence in the future.



