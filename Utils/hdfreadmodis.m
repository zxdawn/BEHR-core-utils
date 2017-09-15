function [ data, fill_value, scale_factor, offset ] = hdfreadmodis( filename, dsetname, varargin )
%HDFREADMODIS Wrapper around HDFREAD that handles fills, scale, and offset for MODIS files
%   DATA = HDFREADMODIS( FILENAME, DSETNAME ) - Reads in the dataset from
%   FILENAME at internal path DSETNAME. Uses HDFREAD internally as
%   HDFREAD(FILENAME, DSETNAME), so the DSETNAME must be formatted such
%   that HDFREAD understands it. Fill values defined by the _FillValue
%   attribute are converted to NaNs and the scale factor and offset are
%   applied as DATA = scale_factor * (hdf_values - add_offset).
%
%   DATA = HDFREADMODIS( ___, 'fill_crit', FILL_CRIT ) allows you to
%   specify the relative difference allowed between a value in the data set
%   and the fill value for the value to be considered a fill value, i.e.
%
%       abs((val - fill_val)/fill_val) < fill_crit 
%
%   must be met for VAL to be considered a fill value. Defaults to 0.001.
%
%   [ DATASET, FILL_VALUE, SCALE_FACTOR, OFFSET ] = H5READOMI( __ ) returns
%   the other attributes as well as the dataset.

E = JLLErrors;
p = inputParser;
p.addParameter('fill_crit',0.001);
p.parse(varargin{:});
pout = p.Results;

fill_crit = pout.fill_crit;
if ~isnumeric(fill_crit) || ~isscalar(fill_crit) || fill_crit <= 0
    E.badinput('FILL_CRIT must be a positive, scalar number')
end

data = double(hdfread(filename, dsetname));
fill_value = double(hdfreadatt(filename, dsetname, '_FillValue'));
try
    scale_factor = double(hdfreadatt(filename, dsetname, 'scale_factor'));
catch err
    if strcmp(err.identifier, 'hdfreadatt:hdf_attribute')
        warning('No scale_factor attribute found for %s, setting to 1', dsetname);
        scale_factor = 1;
    else
        rethrow(err)
    end
end

try
    offset = double(hdfreadatt(filename, dsetname, 'add_offset'));
catch err
    if strcmp(err.identifier, 'hdfreadatt:hdf_attribute')
        warning('No add_offset attribute found for %s, setting to 1', dsetname);
        offset = 1;
    else
        rethrow(err);
    end
end
fills = abs((data - fill_value) ./ fill_value) < fill_crit;
data(fills) = nan;

% This treatment is given by the MODIS MOD06 theoretical basis document, p.
% 87 (https://modis-atmos.gsfc.nasa.gov/_docs/C6MOD06OPUserGuide.pdf) and
% the metadata for MCD43C1 v006
% (https://ladsweb.modaps.eosdis.nasa.gov/api/v1/filespec/product=MCD43C1&collection=6,
% search BRDF_Albedo_Parameter).
data = (data - offset) * scale_factor;

end

