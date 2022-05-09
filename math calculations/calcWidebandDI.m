function [widebandDI] = calcWidebandDI(narrowbandDF)
%Calculates Wideband Directivity Index [dB]
% Input:
%   narrowbandDF - array of narrowband directivity factors
% Output:
%   widebandDI - wideband directivity index [dB]

widebandDI = 10*log10( 1 / mean(1./narrowbandDF) );

end

