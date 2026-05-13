classdef energy_consumption

    methods (Static)

        % ---- Operation counts ----------------------------------------

        function [NM, NA] = fft_search_counts(L, Nfft, blind)
            if blind
                NM_form = 12 * L;
                NA_form = 9  * L;
            else
                NM_form = 4 * L;
                NA_form = 3 * L;
            end
            NM_fft    = 2 * L * log2(Nfft / L) + 2 * Nfft * log2(L);
            NA_fft    = 3 * L * log2(Nfft / L) + 3 * Nfft * log2(L);
            NM_search = 2 * Nfft;
            NA_search = 2 * Nfft;
            NM_interp = 5;
            NA_interp = 4;
            NM = NM_form + NM_fft + NM_search + NM_interp;
            NA = NA_form + NA_fft + NA_search + NA_interp;
        end

        function [NM, NA] = differential_kay_counts(L)
            NM = 7 * L + 1;
            NA = 4 * L - 2;
        end

        function [NM, NA] = viterbi_counts(N)
            NM = 16 * N + 1;            % 12N + 4N + 1
            NA = 14 * N - 1;            % 9N + 3N + 2(N-1) + 1
        end

        function [NM, NA] = pilots_only_counts()
            NM = 1;
            NA = 1;
        end

        % ---- Energy per bit [fJ] -------------------------------------
        %
        % norm_len : number of symbols over which counts are amortised
        %            (SubframeLen for data-aided FR; observation window D
        %             for blind FR; BlockLen for CR)
        % EAdd_fJ, EMult_fJ : per-bit energy coefficients passed to
        %            energy.receiver  (E_A = EAdd*n, E_M = EMult*n^2)
        % M            : modulation order
        % Oversampling : f/Rs  (use 1 for symbol-rate stages)
        % n            : word length in bits

        function E = fft_search_energy(L, Nfft, blind, norm_len, ...
                                       EAdd_fJ, EMult_fJ, M, Oversampling, n)
            [NM, NA] = energy_consumption.fft_search_counts(L, Nfft, blind);
            NM = NM / norm_len;
            NA = NA / norm_len;
            E  = energy.receiver(NA, NM, EAdd_fJ, EMult_fJ, M, Oversampling, n);
        end

        function E = differential_kay_energy(L, norm_len, ...
                                             EAdd_fJ, EMult_fJ, M, Oversampling, n)
            [NM, NA] = energy_consumption.differential_kay_counts(L);
            NM = NM / norm_len;
            NA = NA / norm_len;
            E  = energy.receiver(NA, NM, EAdd_fJ, EMult_fJ, M, Oversampling, n);
        end

        function E = viterbi_energy(N, EAdd_fJ, EMult_fJ, M, Oversampling, n)
            % N is the block length; counts are amortised over N symbols.
            [NM, NA] = energy_consumption.viterbi_counts(N);
            NM = NM / N;
            NA = NA / N;
            E  = energy.receiver(NA, NM, EAdd_fJ, EMult_fJ, M, Oversampling, n);
        end

        function E = pilots_only_energy(BlockLen, EAdd_fJ, EMult_fJ, M, Oversampling, n)
            [NM, NA] = energy_consumption.pilots_only_counts();
            NM = NM / BlockLen;
            NA = NA / BlockLen;
            E  = energy.receiver(NA, NM, EAdd_fJ, EMult_fJ, M, Oversampling, n);
        end

    end

end
