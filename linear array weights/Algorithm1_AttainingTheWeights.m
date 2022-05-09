function [weights, Zopt, DF, WNG, BW, phi, f_min, beta_opt, L_ind_opt] = Algorithm1_AttainingTheWeights(c, freqs, x_coor, desired_bw, beta_res, L_support_option, use_continuous_kaiser, use_trapezoidal_integration)
%Calulates the beamformers weights according to Algorithm 1.
%Inputs:
%   c - speed of propogation [m/s]. Example: speed of sound in air 343 m/s
%   freqs - array of frequencies [Hz] to compute at. Example: 0:10:8e3
%   x_coor - x-coordinates of array [m]. Example: -0.5:0.1:0.5
%   desired_bw - desired beamwidth [deg]. Example: 15
%   beta_res - resolution of Kaiser window shape factor. Example: 1e-3
%   L_support_option - string. Choose what the active sensors are. Options: 'single', 'supports', 'custom'.
%   use_continuous_kaiser - boolean. If to sample the continuous Kaiser window. Otherwise, uses the discrete Kaiser window.
%   use_trapezoidal_integration - boolean. If to use the trapezoidal integration technique.
%Outputs:
%   weights - the beamformers weights
%   Zopt - the beampattern
%   DF - the narrowband directivity factor
%   WNG - the white noise gain
%   BW - the beamwidth
%   phi - azimuth angles that the beampattern was evaluated at
%   f_min - minimum frequency that attains the desired beamwidth
%   beta_opt - the Kaiser window shape factor
%   L_ind_opt - the Kaiser window support

x_coor = x_coor(:); % Make sure it is a column vector.
M = length(x_coor); % number of sensors
N = (M+1)/2; % number of weights (half of M because it is a symmetric array)
theta_look = pi/2; % elevation angle [rad]. Set to pi/2 for XY plain.
phi = linspace(0,pi,500); % azimuth angles [rad]. We use this resolution for plotting the beampattern.
phi_for_bw = linspace(0,pi,3601); % azimuth angles [rad]. We use a higher resolution here to accurately calculate the beamwidth.

pseudoCoherenceMatrix = calcPseudoCoherenceMatrix(x_coor, zeros(size(x_coor)), freqs, c, false); % We calculate it here to save computations.
steeringVector_for_bw = zeros(M, length(phi_for_bw), length(freqs)); % We calculate it here to save computations.
for f_ind = 1:length(freqs)
    f = freqs(f_ind);
    k = f*2*pi/c;
    steeringVector_for_bw(:,:,f_ind) = exp(-1j*k*x_coor*cos(phi_for_bw));
end

weights = zeros(N, length(freqs));
Zopt = zeros(length(phi), length(freqs));
DF = zeros(1, length(freqs));
WNG = zeros(1, length(freqs));
BW = zeros(1, length(freqs));
f_min = Inf;
beta_opt = zeros(1, length(freqs));
L_ind_opt = zeros(1, length(freqs));

beta_min_value = 0;
beta_max_value = 10;
beta = beta_min_value:beta_res:beta_max_value;

if strcmp(L_support_option,'single')
    L = x_coor(end); % 2*L is the support of the entire array
elseif strcmp(L_support_option,'supports')
    L = x_coor(N+1:end);
elseif strcmp(L_support_option,'custom') % L can have the sub-supports of the array
    L = linspace(x_coor(N+1), 1.5*x_coor(end), 101); % support must cover atleast the 3 innermost mics (therefore N+1).
end

beta_found = nan(length(L), length(freqs));
for L_ind = 1:length(L) % L_ind is denoted "i" in the article.
    beta_ind = length(beta);

    L_curr = L(L_ind);
    beta_curr = beta(beta_ind);

    for f_ind = length(freqs):-1:1
        h = kaiser_taps(M, x_coor, use_continuous_kaiser, beta_curr, L_curr, use_trapezoidal_integration);
        attained_bw = calc_BW3dB(steeringVector_for_bw(:,:,f_ind), h, phi_for_bw, false, false);
        if attained_bw <= desired_bw
            beta_found(L_ind, f_ind) = Inf; % To mark that the beamwidth is too narrow.
        end
        while attained_bw > desired_bw && beta_ind > 1
            beta_ind = beta_ind - 1;
            beta_curr = beta(beta_ind);

            h = kaiser_taps(M, x_coor, use_continuous_kaiser, beta_curr, L_curr, use_trapezoidal_integration);
            attained_bw = calc_BW3dB(steeringVector_for_bw(:,:,f_ind), h, phi_for_bw, false, false);

            beta_found(L_ind, f_ind) = beta_curr;
            f_min = freqs(f_ind);
        end % End While
    end % End For over frequency
end % End For over L

% Find which L achieves maximal DF per frequency:
DF_per_L = zeros(size(L));
for f_ind = 1:length(freqs)
    f = freqs(f_ind);
    k = f*2*pi/c;
    for L_ind = 1:length(L)
        L_curr = L(L_ind);
        beta_curr = beta_found(L_ind, f_ind);
        if isnan(beta_curr)
            DF_per_L(L_ind) = nan;
        elseif isinf(beta_curr)
            DF_per_L(L_ind) = -Inf;
        else
            h = kaiser_taps(M, x_coor, use_continuous_kaiser, beta_curr, L_curr, use_trapezoidal_integration);
            DF_per_L(L_ind) = 1/(h'*pseudoCoherenceMatrix(:,:,f_ind)*h);
        end
    end
    if sum(isnan(DF_per_L))==length(L) % if the beawidth was always to large, use the narrowest beamwidth
        L_ind_opt(f_ind) = length(L);
        beta_opt(f_ind) = beta(1);
    elseif sum(isinf(DF_per_L))==length(L) % if the beawidth was always too narrow, use the widest beamwidth
        L_ind_opt(f_ind) = 1;
        beta_opt(f_ind) = beta(end);
    else
        [~, DF_max_ind] = max(DF_per_L);
        L_ind_opt(f_ind) = DF_max_ind;
        beta_opt(f_ind) = beta_found(DF_max_ind, f_ind);
    end

    L_curr = L(L_ind_opt(f_ind));
    beta_curr = beta_opt(f_ind);
    h = kaiser_taps(M, x_coor, use_continuous_kaiser, beta_curr, L_curr, use_trapezoidal_integration);

    BW(f_ind) = calc_BW3dB(steeringVector_for_bw(:,:,f_ind), h, phi_for_bw, false, false);
    weights(:,f_ind) = h(N:end);
    steeringVector = exp(-1j*k*x_coor*cos(phi)*sin(theta_look));
    Zopt(:,f_ind) = h'*steeringVector;
    WNG(f_ind) = calcWNG(h);
    DF(f_ind) = 1/(h'*pseudoCoherenceMatrix(:,:,f_ind)*h);
end

end % End Main Function

function h = kaiser_taps(M, x_coor, use_continuous_kaiser, beta_curr, L_curr, use_trapezoidal_integration)
%Gets the beamforming weight vector.
w_r = zeros(M,1);
x_ind_support = find(abs(x_coor) <= L_curr + eps);

if use_continuous_kaiser
    w_r(x_ind_support) = besseli(0,beta_curr*sqrt(1-(x_coor(x_ind_support).'/L_curr).^2)).'/besseli(0,beta_curr);
else
    w_r(x_ind_support) = kaiser(length(x_ind_support), beta_curr);
end

if use_trapezoidal_integration
    spacing = conv(x_coor,[1; 0; -1]/2,'valid');
    w_r = w_r.*[diff(x_coor(end-1:end)); spacing; diff(x_coor(end-1:end))];
end

h = w_r/sum(w_r); % normalize

end % End Kasier Taps Function