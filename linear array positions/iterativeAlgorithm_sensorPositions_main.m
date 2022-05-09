close all; clear; clc; % clean the work enviornment
addpath(genpath('..\..\code')); % add folder to MATLAB path

%% Parameters
c = 343; % speed of sound [m/s]
isULA = false; % should the interelement spacing be uniform
M = 11; % number of sensors
N = (M-1)/2;
desired_bw = 15; % desired beamwidth [deg]
beta_res = 1e-2; % resolution of Kaiser window shape factor. Example: 1e-3
freqs = 0:10:8e3; % array of frequencies [Hz] to compute at
useFinitePrecision = true; % if the positions should have finite precision.
precision = 3; % amount of decimals after the point. Example: 3 for mm precision
L_support_option = 'supports'; % choose what the active sensors are. Options: 'single', 'supports', 'custom'
use_continuous_kaiser = true; % if to sample the continuous Kaiser window. Otherwise, uses the discrete Kaiser window
use_trapezoidal_integration = true; % if to use the trapezoidal integration technique

% Initial Point
if isULA
    % initialize interelement spacing
    x0 = 0.5*c/freqs(end);
else
    % initialize sensor locations
    x0 = nonIterativeAlgorithm_sensorPositions(c, freqs, M, desired_bw, 0.1e-2, 1.36, 0.034);
    x0 = x0(N+2:end); % keep only the positive sensors. The rest are known because the array is symmetric.
    if length(x0)~=N
        error('the number of sensors in the initial point is not N')
    end
    x0 = [0.0380    0.0790    0.1440    0.2970    0.8260];
end

[found_positions, found_directivity] = iterativeAlgorithm_sensorPositions(x0,c,isULA,N,desired_bw,freqs,useFinitePrecision,precision,beta_res,L_support_option,use_continuous_kaiser,use_trapezoidal_integration);
if useFinitePrecision
    found_positions = round(found_positions,precision);
end

%% Helper Functions
function [found_positions, found_directivity] = iterativeAlgorithm_sensorPositions(x0,c,isULA,N,desired_bw,freqs,useFinitePrecision,precision,beta_res,L_support_option,use_continuous_kaiser,use_trapezoidal_integration)
%Iterative algorithm for finding the sensors' positions

x0 = x0(:); % Make sure it is a column vector.

options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolFun',1e-1,'TolX',1e-3);
[x, fval] = fminsearch(@helperFunction, x0, options);
found_positions = sort(abs(x));
found_directivity = -fval;

    function fval = helperFunction(x)
        %Helper Function - evaluates the DI for sensors positioned at x
        if useFinitePrecision
            x = round(x,precision);
        end
        if isULA
            x_LA = (-x*N:x:x*N)';
        else
            x_LA = sort([-x(end:-1:1); 0; x]);
        end
        [~, ~, directivityFactor_LA, ~, BW_LA, ~, ~, ~, ~] = Algorithm1_AttainingTheWeights_robust(c, freqs, x_LA, desired_bw, beta_res, L_support_option, use_continuous_kaiser, use_trapezoidal_integration);

        fval = -calcWidebandDI(directivityFactor_LA); % minus so that we will maximize DI

        % for faster convergence, check the following conditions. If they are not met, harm the objective
        if ~issorted(BW_LA,'descend') || sum(BW_LA<desired_bw)>0
            fval = fval + 5;
        end

    end

end % End of iterativeAlgorithm_sensorPositions

