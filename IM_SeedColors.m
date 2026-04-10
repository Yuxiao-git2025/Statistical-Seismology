function C=IM_SeedColors(n,method,seed)
%RANDCOLORS Generate an n-by-3 RGB colormap (values in [0,1]).
%   C = RANDCOLORS(n) returns n random RGB colors (uniform in [0,1]).
%   C = RANDCOLORS(n, method) method:
%       "uniform" (default) : rand(n,3)
%       "bright"            : brighter colors, avoid too-dark values
%       "distinct"          : well-separated hues (not purely random)
%   C = RANDCOLORS(n, method, seed) sets RNG seed for reproducibility.
%
%   Example:
%       C = randColors(10,"distinct",1);
%       scatter(1:10, ones(1,10), 80, C, 'filled');
    if nargin < 2 || isempty(method), method = "uniform"; end
    if nargin >= 3 && ~isempty(seed), rng(seed); end
    validateattributes(n, {'numeric'}, {'scalar','integer','positive'});
    method = string(method);
    switch lower(method)
        case "uniform"
            C = rand(n,3);

        case "bright"
            % Sample in HSV: random hue, medium-high saturation & value
            H = rand(n,1);
            S = 0.6 + 0.4*rand(n,1);  % [0.6, 1.0]
            V = 0.75 + 0.25*rand(n,1);% [0.75, 1.0]
            C = hsv2rgb([H S V]);

        case "distinct"
            % Evenly spaced hues + small jitter, then convert HSV->RGB
            H = mod((0:n-1)'/n + 0.02*randn(n,1), 1);
            S = 0.65 + 0.25*rand(n,1);
            V = 0.85 + 0.15*rand(n,1);
            C = hsv2rgb([H S V]);

        otherwise
            error('Unknown method: %s. Use "uniform", "bright", or "distinct".', method);
    end
    C = max(0, min(1, C));
end
