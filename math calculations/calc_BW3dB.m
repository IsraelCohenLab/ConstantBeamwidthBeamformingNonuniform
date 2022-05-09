function [attained_bw] = calc_BW3dB(steeringVector, h_bw, theta, elevation, debug)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
beta_bw = sqrt(0.5); % 3dB bandwidth
dtheta = theta(2)-theta(1);

Zopt = h_bw'*steeringVector; % Beampattern
Zopt_abs = abs(Zopt); % [dB]

% precision for when Zopt_abs~beta_bw
[~,max_ind] = max(Zopt_abs);
epss = 1e-4;%( Zopt_abs(max_ind)-Zopt_abs(max_ind+1) ) / 10;

% find 3dB bandwidth
if elevation
    f1 = 0;
    f2 = find(Zopt_abs <= beta_bw - epss, 1, 'first') - 2; % [ind] (-1 because the first ind is theta=0 and another -1 because we want when the BW is greater than beta_bw)
    attained_bw = 2*(f2-f1)*dtheta*180/pi; % [deg]
else
    f1 = find(Zopt_abs(1:ceil(length(Zopt_abs)/2)) <= beta_bw - epss, 1, 'last')+1; % [ind]
    f2 = find(Zopt_abs(floor(length(Zopt_abs)/2):end) <= beta_bw - epss, 1, 'first')+floor(length(Zopt_abs)/2)-2; % [ind]
    attained_bw = (f2-f1)*dtheta*180/pi; % [deg]
end

if isempty(attained_bw)
    attained_bw = 180;
end

if debug
    figure()
    plot(theta*180/pi,Zopt_abs)
    title('Beam Pattern')
    xlabel('\theta [deg]')
    ylabel('|Beam Pattern|')
    if attained_bw~=180
        xline((f1-1)*dtheta*180/pi); % -1 because indices start at 1
        xline((f2-1)*dtheta*180/pi); % -1 because indices start at 1
    end
end

