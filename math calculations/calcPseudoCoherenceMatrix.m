function pseudoCoherenceMatrix = calcPseudoCoherenceMatrix(r_coor, phi_coor, freqs, c, isPolar)
% calculate Pseudo Coherence Matrix
%   freqs - array of frequencies to calculate matrix [Hz]
%   c - speed of sound [m/s]. For example: 343
%   isPolar - are the coordinates polar

%convert polar to cartesian coordinates
if isPolar
    x_coor = r_coor.*cos(phi_coor);
    y_coor = r_coor.*sin(phi_coor);
else
    x_coor = r_coor;
    y_coor = phi_coor;
end

% the matrix is symmetric with ones along the diagonal
pseudoCoherenceMatrix = ones(length(x_coor), length(x_coor), length(freqs));
for f_ind = 1:length(freqs)
    f = freqs(f_ind);
    for m = 1:length(x_coor)-1
        for n = m+1:length(x_coor)
            pseudoCoherenceMatrix(m,n,f_ind) = sinc(2*f/c*sqrt( (y_coor(m)-y_coor(n))^2 + (x_coor(m)-x_coor(n))^2 ));
            pseudoCoherenceMatrix(n,m,f_ind) = pseudoCoherenceMatrix(m,n,f_ind);
        end
    end
end
