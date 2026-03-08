function Y = add_awgn(X, SNR)
%ADD_AWGN  Add white Gaussian noise to a signal.
%
%   Y = add_awgn(X, SNR)
%
%   Inputs
%     X   - input signal [samples x N_pol]
%     SNR - signal-to-noise ratio [dB]

    Y = awgn(X, SNR, 'measured');
end
