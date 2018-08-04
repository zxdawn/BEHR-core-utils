% Xin 16 May 2018: Only support lno and lnox

% INTEGPR2_lnox Integrate mixing ratios over pressure
%
%   VCD = INTEGPR2_lnox(MIXINGRATIO, PRESSURE, PRESSURESURFACE) Integrates the
%   profile given by MIXINGRATIO (which must be unscaled, i.e. mol/mol or
%   volume/volume) defined at the vertical levels given by PRESSURE (which
%   must be monotonically decreasing and given in hPa). PRESSURESURFACE is
%   the lower level of integration, must be a scalar also given in hPa. The
%   upper limit is taken as the minimum value in PRESSURE.

% Legacy comment block from Ashley Russell:
%..........................................................................
% Integrates vector of mixing ratios (cm-3) above pressureSurface as function of pressure (hPa)
% to get vertical column densities (cm-2). Computes piecewise in layers between two pressures.
% Assumes exponential variation between mixing ratios that are positive or zero (zero is treated as 1.e-30).
% If one or both mixing ratios is negative, assumes constant (average=(f1+f2)/2) mixing ratio in the layer.
% The uncertainties are computed assuming constant mixing ratio in each layer (using exponential
% variation can lead to strange (large) results when f1*p1 is approx f2*p2. This sometimes happens
% when using averaging kernels).
%
% Arguments (vector indices are i = 0,1,2,...n-1):
%
%  mixingRatio(i)   = input volume mixing ratios (no units)
%  pressure(i)      = vector of input pressures from largest to smallest (hPa)
%  pressureSurface  = surface pressure (where integration starts) (hPa)
%  interpPres       = pressures to interpolate the output profile to (hPa) ( added by JLL for BEHR scattering weights )
%
% Restrictions:
%  Pressures must be greater than zero and monotonically decreasing
%
% Equation for exponential variation:
%  For a segment of the vertical column between ambient pressures p1 and p2,
%  at which the trace gas volume mixing ratios are f1 and f2, the vertical
%  column density VCD (assuming number density varies exponentially with p) is
%
%    VCD   =  (f1 * p1  -  f2 * p2) /  [(b + 1) * m * g]
%
%    where:   b = ALOG( f2/f1 ) / ALOG( p2/p1 )
%
%   Also, for interpolation of f between p1 and p2:
%
%   f  =  f1 * (p/p1)^b
%
% Note:  Can get large uncertainties when  f1*p1 = f2*p2
%
%..........................................................................

function [vcd, vcd_lno2] = integPr2_lnox(mixingRatio_lno, mixingRatio_lno2, pressure, pressureSurface, varargin)

E = JLLErrors;
p = inputParser;
p.addOptional('pressureTropopause', [], @(x) isscalar(x) && isnumeric(x) && (isnan(x) || x > 0));
p.addParameter('interp_pres', []);
p.addParameter('fatal_if_nans', false);

p.parse(varargin{:});
pout = p.Results;

pressureTropopause = pout.pressureTropopause;
interpPres = pout.interp_pres;
fatal_if_nans = pout.fatal_if_nans;

if any(pressure<0)
    E.badinput('PRESSURE must be all >= 0')
elseif any(diff(pressure)>0)
    E.badinput('PRESSURE must be monotonically decreasing')
end

if ~isscalar(pressureSurface)
    E.badinput('PRESSURESURFACE must be a scalar')
end

if isempty(pressureTropopause)
   pressureTropopause = min(pressure);
end

if isempty(interpPres)
    if nargout > 1
        E.callError('nargout','Without any interpPres values, p_out and f_out will not be set');
    end
else
    % Make sure the interpolation pressure is within the pressure vector
    % given.
    interpPres = clipmat(interpPres, min(pressure), max(pressure));
end


%   mean molecular mass (kg)  *  g (m/s2)  *  (Pa/hPa)   *   (m2/cm2)
mg_lno  = (30.01/6.02E23)*1E-3    *   9.8      *    1E-2     *     1E4;
mg_lno2  = (28.97/6.02E23)*1E-3    *   9.8      *    1E-2     *     1E4;

fmin = 1E-30;

vcd  = 0;
vcd_lno2 = 0;
f_lno    = mixingRatio_lno;
f_lno2    = mixingRatio_lno2;
p    = pressure;
p0   = max(pressureSurface, min(pressure));
n    = numel(p);

dvcd     = 0;
deltaVcd_lno = zeros(numel(f_lno),1); % Changed to make these vectors on 9/26/2014 JLL
deltaVcd_lno2 = zeros(numel(f_lno2),1); % Changed to make these vectors on 9/26/2014 JLL

% If the surface pressure is above the tropopause pressure, i.e. the lower
% integration limit is above the upper integration limit, return 0 b/c
% unlike abstract integration where reversing the integration limits
% reverses the sign, here, physically, if the surface pressure is above the
% top limit, then there is no column between the surface and the top limit.
% With this added here, we might want to remove the "pCld > pTropo" tests
% in omiAmfAK2.
if pressureSurface <= pressureTropopause
    vcd = 0;
    vcd_lno2 = 0;
    return
end

f0_lno = interpolate_surface_pressure(p,f_lno,p0);
f0_lno2 = interpolate_surface_pressure(p,f_lno2,p0);

p(i0) = p0;
f_lno(i0) = f0_lno;
f_lno2(i0) = f0_lno2;
i_initial = i0;

if ~isnan(pressureTropopause)
    p_end   = max(pressureTropopause, min(pressure));
    f_lno_end = interpolate_surface_pressure(p,f_lno,p_end);
    f_lno2_end = interpolate_surface_pressure(p,f_lno2,p_end);
    p(i0+1) = p_end;
    f_lno(i0+1) = f_lno_end;
    f_lno2(i0+1) = f_lno2_end;
    i_end = i0;
else 
    i_end = n-1;
end


% Integrate................................................................
if isnan(pressureSurface)
    vcd = nan;
    vcd_lno2 = nan;
    return
end

for i = 1:n-1
    deltaVcd_lno(i) = (f_lno(i) + f_lno(i + 1)) .* (p(i) - p(i + 1)) ./ (2 * mg_lno);  %assume const mixing ratio in each layer
    deltaVcd_lno2(i) = (f_lno2(i) + f_lno2(i + 1)) .* (p(i) - p(i + 1)) ./ (2 * mg_lno2);  %assume const mixing ratio in each layer
    b_lno = (log(max(f_lno(i + 1),fmin)) - log(max(f_lno(i),fmin))) ./ (log(p(i + 1)) - log(p(i)));
    b_lno2 = (log(max(f_lno2(i + 1),fmin)) - log(max(f_lno2(i),fmin))) ./ (log(p(i + 1)) - log(p(i)));

    if f_lno(i) >= 0 && f_lno(i + 1)>=0 && abs(b_lno + 1) >= 0.01;
        deltaVcd_lno(i) = (f_lno(i)*p(i) - f_lno(i + 1)*p(i + 1)) ./ (b_lno + 1) ./ mg_lno;    %assume exponential variation in each layer
    end
    if f_lno2(i) >= 0 && f_lno2(i + 1)>=0 && abs(b_lno2 + 1) >= 0.01;
        deltaVcd_lno2(i) = (f_lno2(i)*p(i) - f_lno2(i + 1)*p(i + 1)) ./ (b_lno2 + 1) ./ mg_lno2;    %assume exponential variation in each layer
    end    
end

if any(isnan(deltaVcd_lno(i_initial:i_end)))
    if fatal_if_nans
        E.callError('nans_in_column', 'NaNs detected in partial columns.')
    else
        warning('integPr2:nans_in_column', 'NaNs detected in partial columns. They will not be added into the total column density.')
    end
end

vcd = nansum2(deltaVcd_lno(i_initial:i_end)+deltaVcd_lno2(i_initial:i_end));
vcd_lno2 = nansum2(deltaVcd_lno2(i_initial:i_end));


    function [f0] = interpolate_surface_pressure(p_in,f_in,p0_in)
        pp=find(p_in>=p0_in); if isempty(pp); pp=0; end
        i0 = min(max((max(pp)),1),(n - 1));
        f0 = interp1(log([p_in(i0),p_in(i0 + 1)]), [f_in(i0),f_in(i0 + 1)], log(p0_in),'linear','extrap');   %assume linear variation in surface layer
        
        b = (log(max(f_in(i0 + 1),fmin)) - log(max(f_in(i0),fmin))) ./ (log(p_in(i0 + 1)) - log(p_in(i0)));
        if f_in(i0) >= 0 && f_in(i0 + 1) >= 0 && abs(b + 1) >= 0.01;
            f0 = f_in(i0) * ((p0_in./p_in(i0))^b);                                    %assume exponential variation in surface layer
        end
        
    end
end
