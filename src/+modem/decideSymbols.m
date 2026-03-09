function decidedSymbols = decideSymbols(rxSymbols)
%DECIDESYMBOLS  QPSK hard-decision: nearest ±1 ±1j constellation point.
%
%   decidedSymbols = decideSymbols(rxSymbols)
%
%   Input
%     rxSymbols - [Nsym x NPol] complex received symbols
%
%   Output
%     decidedSymbols - same size, each element snapped to ±1 ±1j

    d_re = sign(real(rxSymbols));
    d_im = sign(imag(rxSymbols));

    d_re(d_re == 0) = 1;
    d_im(d_im == 0) = 1;

    decidedSymbols = complex(d_re, d_im);
end
