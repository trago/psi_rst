% This function generates the synthetic phase used in numerical tests of
% the article: "Regularized self-tuning phase demodulation from a
% phase-shifted interferogram sequence with arbitrary phase shifts".

% Date: June 4 2011
% Author: Orlando M Medina Cazares.

% ----- Inputs ------
% M : Number of rows.
% N : Number of columns.
% A : Amplitud of the phase.

% ----- Return ------
% phase : Synthetic phase generated of M by N pixels.


function phase = makePeaks(M,N,A)
sx = 2*pi/(M-1);
sy = 2*pi/(N-1);
[X Y]   = meshgrid(-pi:sx:pi,-pi:sy:pi); 
phase = -(1 - X/2 + X.^5 + Y.^3) .* exp(-X.^2 - Y.^2 )*A;
end