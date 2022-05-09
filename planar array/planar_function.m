function [weights_PA,Zopt_PA,DF_PA,WNG_PA,BW_xz,BW_yz,PA_x,PA_y,phi,theta] = planar_function(LA_A,LA_B,desired_bw_xz,desired_bw_yz,freqs,f_plot,c,planar_pattern,alpha,do_kronecker)
%Creates a planar array based on linear arrays
%Inputs:
%   LA_A - array of x-coordinates [m] of the planar grid
%   LA_B - array of y-coordinates [m] of the planar grid
%   desired_bw_xz - desired beamwidth [deg] in the XZ plain. Example: 15
%   desired_bw_yz - desired beamwidth [deg] in the YZ plain. Example: 15
%   freqs - array of frequencies [Hz] to compute at. Example: 0:10:8e3
%   f_plot - array of frequencies [Hz] to plot beampatterns at
%   c - speed of propogation [m/s]. Example: speed of sound in air 343 m/s
%   planar_pattern - which sensors of the grid to use. For examples, see getSupportOfPA.m
%   alpha - trade-off parameter
%   do_kronecker - if to use Kronecker product. Otherwise, uses trade-off beamformer
%Outputs:
%   weights_PA - the beamformers weights
%   Zopt_PA - the beampattern
%   DF_PA - the narrowband directivity factor
%   WNG_PA - the white noise gain
%   BW_xz - the beamwidth in the XZ plain
%   BW_yz - the beamwidth in the YZ plain
%   PA_x - the x-coordinates of the planar array's sensors
%   PA_y - the y-coordinates of the planar array's sensors
%   phi - azimuth angles that the beampattern was evaluated at
%   theta - elevation angles that the beampattern was evaluated at

LA_A = LA_A(:); % Make sure it is a column vector.
LA_B = LA_B(:); % Make sure it is a column vector.

%% Parameters
% weight parameters
use_continuous_kaiser = true; % if to sample the continuous Kaiser window. Otherwise, uses the discrete Kaiser window
use_trapezoidal_integration = true; % if to use the trapezoidal integration technique
L_support_option = 'supports'; % choose what the active sensors are. Options: 'single', 'supports', 'custom'
beta_res = 1e-3; % resolution of Kaiser window shape factor. Example: 1e-3
closed_form = false; % if to use the closed form solution for the trade-off beamformer. Otherwise, uses MATLAB quadprog to find the solution

%% Get Linear Arrays Weights
M_A = length(LA_A); % amount of sensors in the linear array A
M_B = length(LA_B); % amount of sensors in the linear array B

weights_LA_A = Algorithm1_AttainingTheWeights_robust(c, freqs, LA_A, desired_bw_xz, beta_res, L_support_option, use_continuous_kaiser, use_trapezoidal_integration);
weights_LA_B = Algorithm1_AttainingTheWeights_robust(c, freqs, LA_B, desired_bw_yz, beta_res, L_support_option, use_continuous_kaiser, use_trapezoidal_integration);

h_A = [weights_LA_A(end:-1:2,:); weights_LA_A];
h_B = [weights_LA_B(end:-1:2,:); weights_LA_B];
h_LAs = [h_A; h_B];

%% Design Planar Array
% S - M_B x M_A boolean support matrix - indicator if there is a sensor at each grid position.
S = getSupportOfPA(M_B, M_A, planar_pattern);

M_PA = sum(S(:)); % amount of sensors in the planar array
[PA_x, PA_y, PA_z] = getCoorPA(S, LA_A, LA_B); % get the coordinates based on the grid.

if ~do_kronecker
    [~, A_support] = transformWeightsToPA(S, M_A, M_B); % transforms matrix: weights calculated for the linear arrays -> weights for the PA

    pseudoCoherenceMatrix = calcPseudoCoherenceMatrix(PA_x, PA_y, freqs, c, false);
    pCoherenceAlpha = (1-alpha) * pseudoCoherenceMatrix + alpha * repmat(eye(M_PA),[1,1,length(freqs)]);
end

%% Get Beampattern
phi = linspace(0,2*pi,361); % azimuth angles [rad]
theta = linspace(0,pi,181); % elevation angle [rad]
theta_for_bw = linspace(0,pi/2,1801); % azimuth angles [rad]. We use a higher resolution here to accurately calculate the beamwidth.

[PHI, THETA] = meshgrid(phi, theta);
PHI = PHI(:).';
THETA = THETA(:).';

weights_PA = zeros(M_PA, length(freqs));
Zopt_PA = zeros(length(theta), length(phi), length(freqs));
DF_PA = zeros(length(freqs),1);
WNG_PA = zeros(length(freqs),1);
BW_xz = zeros(length(freqs),1); % BW in XZ Plain
BW_yz = zeros(length(freqs),1); % BW in YZ Plain
for f_ind = 1:length(freqs)
    if do_kronecker
        % Kronecker product beamformer
        wPA = kron(h_A(:,f_ind).', h_B(:,f_ind));
        wPA = wPA(:);
    else
        % trade-off beamformer
        if closed_form
            pCoherenceInv_A = pCoherenceAlpha(:,:,f_ind) \ A_support';
            wPA = pCoherenceInv_A * ( (A_support*pCoherenceInv_A) \ h_LAs(:,f_ind) );
        else
            % solves the problem well when the matrice has a high condition number
            wPA = quadprog(pCoherenceAlpha(:,:,f_ind), zeros(M_PA, 1), [], [], A_support, h_LAs(:,f_ind));
        end
    end    

    weights_PA(:,f_ind) = wPA/sum(wPA); % it is already normalized...
    
    f = freqs(f_ind);
    k = f*2*pi/c;
    DF_PA(f_ind) = calcDirectivityFactor(weights_PA(:,f_ind), PA_x, PA_y, f, false);
    WNG_PA(f_ind) = calcWNG(weights_PA(:,f_ind));

    %calculate BW for the slice phi=0 (XZ plain) and phi=pi/2 (YZ plain)
    phi_xz = 0; % [rad]
    steeringVector = exp(-1j*k*(PA_x*(cos(phi_xz).*sin(theta_for_bw)) + PA_y*(sin(phi_xz).*sin(theta_for_bw)) + PA_z*cos(theta_for_bw)));
    BW_xz(f_ind) = calc_BW3dB(steeringVector, weights_PA(:,f_ind), theta_for_bw, true, false);
    
    phi_yz = pi/2; % [rad]
    steeringVector = exp(-1j*k*(PA_x*(cos(phi_yz).*sin(theta_for_bw)) + PA_y*(sin(phi_yz).*sin(theta_for_bw)) + PA_z*cos(theta_for_bw)));
    BW_yz(f_ind) = calc_BW3dB(steeringVector, weights_PA(:,f_ind), theta_for_bw, true, false);
    
    if ismember(f, f_plot)
        steeringVector = exp(-1j*k*(PA_x*(cos(PHI).*sin(THETA)) + PA_y*(sin(PHI).*sin(THETA)) + PA_z*cos(THETA)));
        Z = weights_PA(:,f_ind)'*steeringVector;
        Z = reshape(Z,length(theta),length(phi));
        Zopt_PA(:,:,f_ind) = Z;
    end
end
