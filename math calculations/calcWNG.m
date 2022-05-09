function WNG = calcWNG(h)
% Calculates the White Noise Gain (WGN). Assumes distortionless constraint.

WNG = 1/(h'*h);
