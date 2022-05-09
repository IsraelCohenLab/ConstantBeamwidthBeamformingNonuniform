function [weights, Zopt, DF, WNG, BW, phi, f_min, beta_opt, L_ind_opt] = Algorithm1_AttainingTheWeights_robust(c, freqs, x_coor, desired_bw, beta_res, L_support_option, use_continuous_kaiser, use_trapezoidal_integration)
%Calulates the beamformers weights according to Algorithm 1. Robust version (works with low resolution of beta).
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

for L_ind = 1:length(L) % L_ind is denoted "i" in the article.
    % initialize with widest beampattern possible
    beta_ind = length(beta);
    h = kaiser_taps(M, x_coor, use_continuous_kaiser, beta(beta_ind), L(L_ind), use_trapezoidal_integration);

    f_ind = length(freqs);
    while f_ind >= 1
        f = freqs(f_ind);
        k = f*2*pi/c;
        attained_bw = calc_BW3dB(steeringVector_for_bw(:,:,f_ind), h, phi_for_bw, false, false);
        % as frequency decreases, the bandwidth increases.
        if attained_bw > desired_bw + 0.05 % the test is to assure that the BW isn't too narrow
            if beta_ind <= 1
                % we can no longer reduce the BW, so we save the best we are able to achieve.
                DF_curr = 1/(h'*pseudoCoherenceMatrix(:,:,f_ind)*h);
                if ( DF_curr >= DF(f_ind) - 10*eps || attained_bw < BW(f_ind) )% -eps because for f=0 there is inaccuracy, and I want the last L to be chosen
                    BW(f_ind) = attained_bw;
                    weights(:,f_ind) = h(N:end);
                    steeringVector = exp(-1j*k*x_coor*cos(phi)*sin(theta_look));
                    Zopt(:,f_ind) = h'*steeringVector;
                    WNG(f_ind) = calcWNG(h);
                    DF(f_ind) = DF_curr;
                    beta_opt(f_ind) = beta(beta_ind);
                    L_ind_opt(f_ind) = L_ind;
                end
                f_ind = f_ind - 1;
            else
                beta_ind = beta_ind - 1;
                h = kaiser_taps(M, x_coor, use_continuous_kaiser, beta(beta_ind), L(L_ind), use_trapezoidal_integration);
            end
        elseif attained_bw >= desired_bw % the test is to assure that the BW isn't too narrow
            DF_curr = 1/(h'*pseudoCoherenceMatrix(:,:,f_ind)*h);
            if ( DF_curr > DF(f_ind) || attained_bw < BW(f_ind) ) % the second test is to assure that the BW isn't too narrow because of a previous L
                f_min = freqs(f_ind);
                BW(f_ind) = attained_bw;
                weights(:,f_ind) = h(N:end);
                steeringVector = exp(-1j*k*x_coor*cos(phi)*sin(theta_look));
                Zopt(:,f_ind) = h'*steeringVector;
                WNG(f_ind) = calcWNG(h);
                DF(f_ind) = DF_curr;
                beta_opt(f_ind) = beta(beta_ind);
                L_ind_opt(f_ind) = L_ind;
            end
            if beta_ind <= 1
                f_ind = f_ind - 1;
            else
                beta_ind = beta_ind - 1;
                h = kaiser_taps(M, x_coor, use_continuous_kaiser, beta(beta_ind), L(L_ind), use_trapezoidal_integration);
            end
        else
            % the BW is narrow
            if beta_ind == length(beta) % if the widest BW is too narrow (for example at high freqs), then we will save the widest BW
                % we can no longer reduce the BW, so we save the widest we are able to achieve.
                if attained_bw > BW(f_ind)
                    DF_curr = 1/(h'*pseudoCoherenceMatrix(:,:,f_ind)*h);
                    BW(f_ind) = attained_bw;
                    weights(:,f_ind) = h(N:end);
                    steeringVector = exp(-1j*k*x_coor*cos(phi)*sin(theta_look));
                    Zopt(:,f_ind) = h'*steeringVector;
                    WNG(f_ind) = calcWNG(h);
                    DF(f_ind) = DF_curr;
                    beta_opt(f_ind) = beta(beta_ind);
                    L_ind_opt(f_ind) = L_ind;
                end
            else
                beta_ind = beta_ind + 1; % the beta is too small (yet previous beta might have been good - and also might be good for the next freq)

                h = kaiser_taps(M, x_coor, use_continuous_kaiser, beta(beta_ind), L(L_ind), use_trapezoidal_integration);
            end
            f_ind = f_ind - 1;
        end
    end
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