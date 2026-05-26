function y = apply_adaptive_eq(In, SpS, opts)
%APPLY_ADAPTIVE_EQ  Run adaptive_eq.equalize from a settings struct.
%
%   y = apply_adaptive_eq(In, SpS, opts)
%
%   Thin wrapper that fills in adaptive_eq.equalize's many optional
%   arguments from an opts struct so the combined eq+clk blocks don't
%   need to hard-code values.  Fields and defaults:
%
%       NTaps          - FIR length        (required)
%       Mu             - step size         (required)
%       SingleSpike    - true
%       N1             - 1
%       NOut           - 0
%       SignOnly       - false
%       PLanes         - 1
%       Mode           - 0   (0 = CMA, 1 = pilot-aided)
%       Pilots         - []  (NBlocks x 2)
%       BlockLen       - PLanes  (when left empty/missing)
%       SubframeBlocks - 0

    o = defaultOpts();
    if nargin >= 3 && ~isempty(opts)
        f = fieldnames(opts);
        for i = 1:numel(f)
            o.(f{i}) = opts.(f{i});
        end
    end
    if isempty(o.BlockLen)
        o.BlockLen = o.PLanes;
    end
    if isempty(o.NTaps) || isempty(o.Mu)
        error('eq_clk:apply_adaptive_eq:missingField', ...
              'opts.NTaps and opts.Mu are required.');
    end

    y = adaptive_eq.equalize(In, SpS, o.NTaps, o.Mu, o.SingleSpike, ...
        o.N1, o.NOut, o.SignOnly, o.PLanes, o.Mode, o.Pilots, ...
        o.BlockLen, o.SubframeBlocks);
end


function o = defaultOpts()
    o.NTaps          = [];
    o.Mu             = [];
    o.SingleSpike    = true;
    o.N1             = 1;
    o.NOut           = 0;
    o.SignOnly       = false;
    o.PLanes         = 1;
    o.Mode           = 0;
    o.Pilots         = [];
    o.BlockLen       = [];
    o.SubframeBlocks = 0;
end
