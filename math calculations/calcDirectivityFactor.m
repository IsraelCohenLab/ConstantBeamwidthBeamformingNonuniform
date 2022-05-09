function directivityFactor = calcDirectivityFactor(h, r_coor, phi_coor, f, isPolar)
% Calculates the Directivity Factor. Assumes distortionless constraint.
%   h - beamformer weights
%   f - frequency [Hz]
%   isPolar - are the coordinates polar

c = 343; % speed of sound [m/s]

%convert polar to cartesian coordinates
if isPolar
    x_coor = r_coor.*cos(phi_coor);
    y_coor = r_coor.*sin(phi_coor);
else
    x_coor = r_coor;
    y_coor = phi_coor;
end

% calculate Pseudo Coherence Matrix (the matrix is symmetric with ones along the diagonal)
gamma = ones(length(x_coor));
for m = 1:length(x_coor)-1
    for n = m+1:length(x_coor)
        gamma(m,n) = sinc(2*f/c*sqrt( (y_coor(m)-y_coor(n))^2 + (x_coor(m)-x_coor(n))^2 ));
        gamma(n,m) = gamma(m,n);
    end
end

directivityFactor = 1/(h'*gamma*h);
