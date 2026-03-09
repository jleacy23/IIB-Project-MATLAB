function bits = symbolsToBits(symbols)
%SYMBOLSTOBITS  Demap QPSK symbols to bits.
%
%   bits = symbolsToBits(symbols)
%
%   Input
%     symbols - [Nsym x NPol] complex QPSK symbols.
%               NPol = 2 : CPON DP-QPSK interleaved bit ordering
%               NPol = 1 : simple I, Q ordering (2 bits per symbol)
%
%   Output
%     bits    - column vector of 0/1
%               NPol = 2 : [Nsym*4 x 1]  (CPON: c(4i)=XI, c(4i+1)=YI,
%                                                 c(4i+2)=XQ, c(4i+3)=YQ)
%               NPol = 1 : [Nsym*2 x 1]  (I then Q per symbol)

    [Nsym, NPol] = size(symbols);

    if NPol >= 2
        % Dual-pol CPON interleaving
        bXI = double(real(symbols(:,1)) > 0);
        bYI = double(real(symbols(:,2)) > 0);
        bXQ = double(imag(symbols(:,1)) > 0);
        bYQ = double(imag(symbols(:,2)) > 0);

        bits = zeros(Nsym * 4, 1);
        bits(1:4:end) = bXI;
        bits(2:4:end) = bYI;
        bits(3:4:end) = bXQ;
        bits(4:4:end) = bYQ;
    else
        % Single-pol: 2 bits per symbol (I, Q)
        bI = double(real(symbols) > 0);
        bQ = double(imag(symbols) > 0);

        bits = zeros(Nsym * 2, 1);
        bits(1:2:end) = bI;
        bits(2:2:end) = bQ;
    end
end
