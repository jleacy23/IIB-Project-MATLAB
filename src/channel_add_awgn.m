function Y = channel_add_awgn(X, SNR)
%CHANNEL_ADD_AWGN  Add white Gaussian noise to a signal.
%
%   Y = channel_add_awgn(X, SNR)
%
%   Inputs
%     X   - input signal [samples x N_pol]
%     SNR - signal-to-noise ratio [dB]

    Y = awgn(X, SNR, 'measured');
end
