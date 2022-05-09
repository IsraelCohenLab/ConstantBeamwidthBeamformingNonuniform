function [x_coor, f_mins, fmin_with_beta0] = nonIterativeAlgorithm_sensorPositions(c, freqs, M, desired_bw, p_res, beta_min, d_init)
%Non-iterative algorithm for finding the sensors' positions
%Inputs:
%   c - speed of propogation [m/s]. Example: speed of sound in air 343 m/s
%   freqs - array of frequencies [Hz] to compute at. Example: 0:10:8e3
%   M - number of sensors. Example: 11
%   desired_bw - desired beamwidth [deg]. Example: 15
%   p_res - resolution of positions [m]. Example: 1e-2
%   beta_min - minimal beta used to create the narrowest beampatterns. Lower beta allows more narrow (therefore lower freqs) but higher sidelobes
%   d_init - the interelement spacing [m] of the initial 5 sensor ULA
%Outputs:
%   x_coor - the sensors' positions
%   f_mins - the minimum frequency that can maintain the desired beamwidth with beta_min (per the amount of sensors used)
%   fmin_with_beta0 - minimum frequency that attains the desired beamwidth

phi_for_bw = linspace(0,pi,3601); % azimuth angles [rad]
beta_max = 10;
N = (M-1)/2;
f_min_inds = zeros(1,N);

% Init m=5 array
m = 5;
n = (m-1)/2;
x_coor = ((-n:n)*d_init).';

f = freqs(end);
k = f*2*pi/c;

w_r = besseli(0,beta_min*sqrt(1-(x_coor.'/x_coor(end)).^2)).'/besseli(0,beta_min);
spacing = conv(x_coor,[1; 0; -1]/2,'valid');
w_r = w_r.*[diff(x_coor(end-1:end)); spacing; diff(x_coor(end-1:end))];
h = w_r/sum(w_r); % normalize

for f_ind = length(freqs):-1:1
    f = freqs(f_ind);
    k = f*2*pi/c;
    steeringVector_for_bw = exp(-1j*k*x_coor*cos(phi_for_bw));
    attained_bw = calc_BW3dB(steeringVector_for_bw, h, phi_for_bw, false, false);
    if attained_bw>=desired_bw
        break
    end
end
f_min_inds(1) = length(freqs);
f_min_inds((m-1)/2) = f_ind;

for n = (m-1)/2+1:N
    x_coor = [x_coor(1)-p_res; x_coor; x_coor(end)+p_res];
    while true
        % narrowest window -  widest BW - too be able to support fmin(n-1)
        w_r = besseli(0,beta_max*sqrt(1-(x_coor.'/x_coor(end)).^2)).'/besseli(0,beta_max);
        spacing = conv(x_coor,[1; 0; -1]/2,'valid');
        w_r = w_r.*[diff(x_coor(end-1:end)); spacing; diff(x_coor(end-1:end))];
        h = w_r/sum(w_r); % normalize

        steeringVector_for_bw = exp(-1j*k*x_coor*cos(phi_for_bw));
        attained_bw = calc_BW3dB(steeringVector_for_bw, h, phi_for_bw, false, false);
        if attained_bw<=desired_bw + 0.05
            % calculate fmin
            w_r = besseli(0,beta_min*sqrt(1-(x_coor.'/x_coor(end)).^2)).'/besseli(0,beta_min);
            spacing = conv(x_coor,[1; 0; -1]/2,'valid');
            w_r = w_r.*[diff(x_coor(end-1:end)); spacing; diff(x_coor(end-1:end))];
            h = w_r/sum(w_r); % normalize
            for f_ind = f_min_inds(n-1):-1:1
                f = freqs(f_ind);
                k = f*2*pi/c;
                steeringVector_for_bw = exp(-1j*k*x_coor*cos(phi_for_bw));
                attained_bw = calc_BW3dB(steeringVector_for_bw, h, phi_for_bw, false, false);
                if attained_bw>=desired_bw
                    break
                end
            end
            f_min_inds(n) = f_ind;
            break
        else
            x_coor(1)   = x_coor(1)   - p_res;
            x_coor(end) = x_coor(end) + p_res;
        end
    end
end

f_mins = freqs(f_min_inds);

%% Calculate f_min when beta = 0
% calculate fmin
w_r = besseli(0,0*sqrt(1-(x_coor.'/x_coor(end)).^2)).'/besseli(0,0);
spacing = conv(x_coor,[1; 0; -1]/2,'valid');
w_r = w_r.*[diff(x_coor(end-1:end)); spacing; diff(x_coor(end-1:end))];
h = w_r/sum(w_r); % normalize
for f_ind = f_min_inds(end):-1:1
    f = freqs(f_ind);
    k = f*2*pi/c;
    steeringVector_for_bw = exp(-1j*k*x_coor*cos(phi_for_bw));
    attained_bw = calc_BW3dB(steeringVector_for_bw, h, phi_for_bw, false, false);
    if attained_bw==desired_bw
        fmin_with_beta0 = freqs(f_ind);
        break
    elseif attained_bw > desired_bw % this freq is already too wide, yet the previous freq was narrow enough
        fmin_with_beta0 = freqs(f_ind+1);
        break
    end
end


