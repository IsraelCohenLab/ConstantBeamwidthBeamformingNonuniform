function [Ainv, A_support] = transformWeightsToPA(S, M_A, M_B)
% transforms the weights calculated for linear arrays to attain weights for a Planar Array
% M_A - amount of sensors in the linear array A. The beamwidth of the PA in the XZ plain will match this array's.
% M_B - amount of sensors in the linear array B. The beamwidth of the PA in the YZ plain will match this array's.
% S - M_B x M_A boolean support matrix - indicator if there is a sensor at each grid position.

% Returns: Ainv
% Use: h_PA = Ainv * [h_A; h_B] (LS solution)
% h_A - weights for linear array A.
% h_B - weights for linear array B.
% make sure h_A and h_B are collumn vectors
% Also valid solution: h_PA = Ainv*[h_A; h_B]+(I-Ainv*A_support)*w  - for any w

S_inds = find(S); % sensors linear indices in the PA

A_complete = zeros(M_A + M_B, M_A * M_B);
for m = 1:M_A
    A_complete(m, (m-1)*M_B+1 : (m-1)*M_B+M_B) = 1;
end
for m = 1:M_B
    A_complete(M_A + m, m:M_B:end) = 1;
end
A_support = A_complete(:,S_inds);

Ainv = pinv(A_support);


