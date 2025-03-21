%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to generate dictionary of frequency-dependent RST (greenFRaySpace) for a grid of point sources 
%
% Marco Olivieri & Mirco Pezzoli 
% 22/06/2023
% v 0.1

clear all, close all, clc;

addpath(genpath('lib'));
addpath(genpath('audio'));


% *******************************************************
%% SETUP

params.c = 343;      % Sound speed
params.fs = 8000;    % Sampling frequency
params.Nfft = 256;   % Number of FFT points
params.winlen = 256; % window lenth
params.ovlap = 128;  % overlap
params.fmin = 20;    % minimum frequency for processing
params.fmax = 2000;  % maximum frequency for processing

% Computed from above
params.fbinsep = params.fs/params.Nfft;
fbins = (0:round(params.Nfft/2))*params.fs/params.Nfft;
minbin = round(params.fmin/params.fbinsep + 1);
maxbin = round(params.fmax/params.fbinsep + 1);
kbins = 2*pi*fbins ./ params.c;
mybins = minbin:maxbin;
nbins = maxbin-minbin+1;


%=====================================================
% Define the array - along y axis
nMic = 32;
micAperure = 3.2;
%micPos = [linspace(0, micAperure, nMic).', zeros(nMic,1), zeros(nMic,1)];
micPos = [zeros(nMic,1), linspace(0, micAperure, nMic).' zeros(nMic,1)];
micPos = micPos.';
d = micPos(2,2);
array.pos = micPos;

% Define a grid of points
gridPx = linspace(0,2,10);
gridPy = linspace(0,micAperure,10);
[XX, YY] = meshgrid(gridPx, gridPy);
gridPts = [XX(:), YY(:), ones(length(gridPx)*length(gridPy),1)*micPos(3,1)].' + [micPos(1,1)+0.3; micPos(1,2); 0] ;

%====================================================
% Ray Space Transform parameters
mubar = 0.06;        % m axis sampling interval (0.6)
D = 51;              % length of m axis (W)
nubar = d;        % nu axis sampling interval (4*d)
sigma = 4*d;        % gaussian window standard deviation (4*d)
mu = tan(linspace(0, pi, D));%((0:mubar:(D-1)*mubar)-((D-1)/2*mubar))';
mu = ((0:mubar:(D-1)*mubar)-((D-1)/2*mubar))';          % [D,1] mu axis

nu = (0:nubar:micPos(2,end))';                        % [L,1] nu axis
t = 1:nMic;

figure
scatter3(array.pos(1,:), array.pos(2,:), array.pos(3,:), 'filled')
hold on
scatter3(gridPts(1,:), gridPts(2,:), gridPts(3,:), 'x')
text(gridPts(1,:)+0.01, gridPts(2,:)+0.01, ...
     gridPts(3,:)+0.01, num2cell(1:length(gridPts)));
xlabel('x'), ylabel('y'), zlabel('z')
axis equal, view(2);
title('Setup');

% *******************************************************
%% RST computation

psi = rayspacetransformmatrix(fbins(mybins),params.c,d,length(micPos),mubar,D,nubar,sigma);
I = size(psi, 1);        % Number of Ray space data points


dist = pdist2(micPos.', gridPts.');     % source microphones distances
for ss = 1:length(gridPts)
    %display(ss)
    fprintf("Simulating point %i / %i \n", ss, length(gridPts))

    % Microphone signals (e.g. free-field Green's Function)
    greenF(ss, :,:) = (exp(-1i*2*pi*fbins(mybins)/params.c .*dist(:,ss)) ./ (4*pi*dist(:,ss))).';  

    % RST computation --> [greenFRaySpace] = [mybins, IxD, gridPts]
    greenFRaySpace(:,:,ss) = squeeze(spatialfilter(permute(greenF(ss, :,:), [2,1,3]), psi,false));
end

% *******************************************************
%% VIEW examples

n_point = 45;

figure
for ff=1:length(mybins)
    % reshape RST vector into image
    Z = reshape(greenFRaySpace(ff,:,n_point), I, D);

    % RST magnitude
    subplot(211)
    imagesc(mu, nu, abs(Z))
    axis xy
    colorbar
    xlabel("m"), ylabel("q"), title("|Z|")

    % RST phase
    subplot(212)
    imagesc(mu, nu, angle(Z))
    axis xy
    colorbar
    xlabel("m"), ylabel("q"), title("\angle Z")

    sgtitle("RST - point " + num2str(n_point) + " - freq " + num2str(round(fbins(ff))) + " Hz")

    pause(0.1)

end


