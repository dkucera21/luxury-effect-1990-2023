%% equity_gap_half_life_and_intersection.m
% Computes (1) the calendar year when the Top–Bottom gap would be HALVED
% and (2) the calendar year when the gap would reach ZERO (lines intersect),
% using FE lines from stacked-tail LMEs:
%   y ~ 1 + time_dec*tail + ...   (Bottom is reference; Top is the 'tail' level)
%
% Requires in workspace: lmeN (NDVI), lmeL (LST) as LinearMixedModel objects.

year0 = 1990;                         % baseline used in time_dec
alpha = 0.05; z = norminv(1 - alpha/2);

assert(exist('lmeN','var')==1 && isa(lmeN,'LinearMixedModel'), 'lmeN not found.');
assert(exist('lmeL','var')==1 && isa(lmeL,'LinearMixedModel'), 'lmeL not found.');

% Build tables for both metrics
T_half = calc_event(lmeN, 'NDVI', year0, z, 'half');   % gap half-life
T_half = [T_half; calc_event(lmeL, 'LST', year0, z, 'half')];

T_int  = calc_event(lmeN, 'NDVI', year0, z, 'zero');   % intersection year
T_int  = [T_int;  calc_event(lmeL, 'LST', year0, z, 'zero')];

disp('— Gap HALF-LIFE (Top–Bottom FE gap reduced to 50%) —');
disp(T_half);
for i = 1:height(T_half)
    feasible = T_half.Direction{i};    % 'closing' / 'widening or parallel'
    yr = T_half.Year(i);  lo = T_half.CI_L(i); hi = T_half.CI_U(i);
    tag = 'outside 1990–2023';
    if isfinite(yr) && yr >= 1990 && yr <= 2023, tag = 'within 1990–2023'; end
    fprintf('%s: half-life ≈ %0.1f (95%% CI %0.1f–%0.1f)  [%s; %s]\n', ...
        T_half.Metric(i), yr, lo, hi, tag, feasible);
end

disp(' ');
disp('— Intersection year (Top vs Bottom FE lines cross; gap → 0) —');
disp(T_int);
for i = 1:height(T_int)
    feasible = T_int.Direction{i};
    yr = T_int.Year(i);  lo = T_int.CI_L(i); hi = T_int.CI_U(i);
    tag = 'outside 1990–2023';
    if isfinite(yr) && yr >= 1990 && yr <= 2023, tag = 'within 1990–2023'; end
    fprintf('%s: intersection ≈ %0.1f (95%% CI %0.1f–%0.1f)  [%s; %s]\n', ...
        T_int.Metric(i), yr, lo, hi, tag, feasible);
end

%% ------------ helpers ------------
function out = calc_event(lme, metric, year0, z, mode)
% mode = 'half' (gap reduced to 50%) or 'zero' (gap = 0)
    CF    = lme.Coefficients;
    cn    = string(CF.Name);                 % term labels
    Sigma = lme.CoefficientCovariance;       % FE covariance

    % Identify the tail (Top vs Bottom) FE and its interaction with time_dec
    iTail = find(contains(cn,'tail') & ~contains(cn,'time_dec'), 1);
    iInt  = find(contains(cn,'time_dec') & contains(cn,'tail'), 1);
    if isempty(iTail) || isempty(iInt)
        error('%s: could not locate tail main effect and/or time_dec:tail interaction.', metric);
    end

    % FE gap model: G(t_dec) = a + b * t_dec
    a = CF.Estimate(iTail);     % baseline FE gap (Top–Bottom) at year0
    b = CF.Estimate(iInt);      % FE change in the gap per decade (Top vs Bottom)

    % By default assume “closing” if a*b < 0, “widening or parallel” otherwise
    if isfinite(a) && isfinite(b) && (a*b < 0)
        direction = "closing";
    elseif abs(b) < 1e-12
        direction = "parallel (no change)";
    else
        direction = "widening or reversing";
    end

    % Solve for t_dec
    switch lower(mode)
        case 'half'         % G(t_half) = 0.5 * a  =>  t = -a/(2b)
            if abs(b) < 1e-12
                t_dec = sign(-a) * inf;
                year  = NaN; ci = [NaN NaN];
            else
                t_dec = -a / (2*b);
                year  = year0 + 10*t_dec;

                % Delta-method SE for t = -a/(2b)
                Va = Sigma(iTail,iTail);
                Vb = Sigma(iInt,iInt);
                Cab= Sigma(iTail,iInt);
                dt_da = -1/(2*b);
                dt_db =  a/(2*b^2);
                Var_t = dt_da^2 * Va + dt_db^2 * Vb + 2*dt_da*dt_db*Cab;
                se_t  = sqrt(max(Var_t,0));
                ci_t  = t_dec + z*[-1 1]*se_t;
                ci    = year0 + 10*ci_t;
            end

        case 'zero'         % G(t_zero) = 0  =>  t = -a/b
            if abs(b) < 1e-12
                t_dec = sign(-a) * inf;
                year  = NaN; ci = [NaN NaN];
            else
                t_dec = -a / b;
                year  = year0 + 10*t_dec;

                % Delta-method SE for t = -a/b
                Va = Sigma(iTail,iTail);
                Vb = Sigma(iInt,iInt);
                Cab= Sigma(iTail,iInt);
                dt_da = -1/b;
                dt_db =  a/(b^2);
                Var_t = dt_da^2 * Va + dt_db^2 * Vb + 2*dt_da*dt_db*Cab;
                se_t  = sqrt(max(Var_t,0));
                ci_t  = t_dec + z*[-1 1]*se_t;
                ci    = year0 + 10*ci_t;
            end

        otherwise
            error('mode must be ''half'' or ''zero''.');
    end

    out = table( ...
        string(metric), string(mode), year, ci(1), ci(2), ...
        a, b, string(direction), ...
        'VariableNames', {'Metric','Mode','Year','CI_L','CI_U','Gap_FE_a','GapSlope_FE_b','Direction'});
end
