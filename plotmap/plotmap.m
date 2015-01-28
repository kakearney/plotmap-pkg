function h = plotmap(varargin)
%PLOTMAP Projects data to a map axis
%
% h = plotmap(param1, val1, param2, val2, ...)
%
% This function creates a map axis and projects the input data to this
% axis, preprocessing the data to determine the proper map and frame
% limits, adding coastlines, and eliminating many of the unecessary data
% restrictions of certain basic mapping functions.  It does not include the
% full functionality of either axesm or the various map plotting functions,
% but eliminates the hassle of creating a simple map, especially if data
% crosses the dateline.
%
% Input variables: 
%
%   latlim:     axis latitude boundaries ([south north])  [-90 90]
%
%   lonlim:     axis longitude boundaries ([west east]).  These can be on
%               either a -180-180 or 0-360 scale. [-180 180] 
%
%   projection: map projection ['mercator']
%
%   line:       n x 2 cell array of data to be plotted via plotm.  Columns
%               represent latitude and longitude, respectively, and each
%               row is plotted as a separate line object.
%
%   contour:    n x 3 cell array of data to be plotted using contourm.
%               Columns represent latitude, longitude, and data,
%               respectively, and each row is plotted as a separate hggroup
%               object.  Currently no flexibility in number of contours.
%
%   pcolor:     n x 3 cell array of data to be plotted using pcolorm.
%               Columns represent latitude, longitude, and data,
%               respectively.  Latitude and longitude can be either vectors
%               or matrices, and data can be any type (unlike restriction
%               to double in pcolorm).  Each row is plotted as a separate
%               surface object.
%
%   patch:      n x 3 cell array of data to be plotted using patchm.
%               Columns represent latitude, longitude, and cdata,
%               respectively, and each row is plotted as a separate patch
%               object.
%
%   patchfv:    n x 3 cell array of data to be plotted using patchm with
%               after first being converted via poly2fv (preserves holes
%               when plotting).  Columns represent latitude, longitude, and
%               cdata, respectively, and each row is plotted as a separate
%               patch object.  
%
%   trim:       If true, calculate longitude and latitude limits based on
%               the extent of data to be plotted and the specified buffer.
%               [true]
%
%   gridint:    1 x 2 vector specifying spacing of meridians and parallels,
%               respectively [30 15]
%
%   coast:      string specifying the GSHHS dataset to use for coastlines:
%               crude (c), low (l), intermediate (i), high(h), or full (f),
%               or use Matlab's world coastlines from coast.mat (m).
%               ['m']
%
%   buffer:     Buffer distance used to calculate frame around data if
%               'trim' is true, in degrees.  With no buffer, points that
%               fall on the boundary will be eliminated when plotting line
%               objects. [1]
%
%   fillcoast:  If true, land is plotted as patch objects.  If false,
%               coastlines are plotted as line objects. [false]
%
%   axis:       handle to existing axis.  If included, map will be added to
%               this axis rather than creating a new figure and axis. 

% Copyright 2008 Kelly Kearney

%-------------------------
% Parse input
%-------------------------

if nargin < 1 || mod(nargin,2)
    error('Input must be passed as parameter-value pairs');
end

In.latlim = [-90 90];
In.lonlim = [-180 180];
In.projection = 'mercator';
In.line = cell(0,2);
In.contour = cell(0,3);
In.pcolor = cell(0,3);
In.patch = cell(0,3);
In.patchfv = cell(0,3);
In.trim = true;
In.gridint = [30 15];
In.coast = 'm';
In.buffer = 1;
In.fillcoast = false;
In.axis = NaN;
In.pv = cell(0);

In = parse_pv_pairs(In, varargin);

% Calculate new coordinates for patchfv objects

if ~isempty(In.patchfv)
    temp = cell(size(In.patchfv));
    temp(:,3) = In.patchfv(:,3);
    for ip = 1:size(In.patchfv, 1)
        [f,v] = poly2fv(In.patchfv{ip,2}, In.patchfv{ip,1});
        vlat = v(:,2);
        vlon = v(:,1);
        temp{ip,1} = vlat(f');
        temp{ip,2} = vlon(f');
    end
    In.patch = [In.patch; temp];
end

if ~isempty(In.patchfv)
    nfv = size(temp,1);
    isfv = false(size(In.patch,1),1);
    isfv(end-(nfv-1):end) = true;
else
    isfv = false(size(In.patch,1),1);
end

%-------------------------
% Calculate map limits
%-------------------------

if In.trim

    % Consider each pcolor grid, each line segment, each individual patch,
    % and each contour grid as a single data object
    
    if ~isempty(In.line)
        latline = cells2vec(In.line(:,1));
        lonline = cells2vec(In.line(:,2));
        [latline, lonline] = polysplit(latline, lonline);
    else
        latline = cell(0);
        lonline = cell(0);
    end
    
    if ~isempty(In.patch)
        latpatch = cells2vec(In.patch(:,1));
        lonpatch = cells2vec(In.patch(:,2));
        [latpatch, lonpatch] = polysplit(latpatch, lonpatch);
    else
        latpatch = cell(0);
        lonpatch = cell(0);
    end

    lon = cat(1, lonline, lonpatch, In.pcolor(:,2), In.contour(:,2));
    lon = cellfun(@changecoords, lon, 'uni', 0);
    
    lat = cat(1, latline, latpatch, In.pcolor(:,1), In.contour(:,1));
    
    ndata = length(lon);
    
    % Now check what the eastern and western bounds of each data object
    % are
    
    west = nan(ndata,1);
    east = nan(ndata,1);
    
    for id = 1:ndata
        
        allinhemis = all(lon{id}(:) >= 0) || all(lon{id}(:) <= 0);
        
        bounds = minmax(lon{id});
        
        if allinhemis
            west(id) = bounds(1);
            east(id) = bounds(2);
        else
            
            iswest = lon{id} <= 0;
            iseast = lon{id} >= 0;
            
            dontknow = isequal(unique(lon{id}), bounds');
            crossesdate = all(lon{id}(iswest) <= bounds(1)) & ...
                          all(lon{id}(iseast) >= bounds(2));
            crossesprime = all(lon{id}(iswest) >= bounds(1)) & ...
                           all(lon{id}(iseast) <= bounds(2));
            
            
            if dontknow
                widthoverprime = diff(bounds);
                widthoverdate = (bounds(1)+180) + (180-bounds(2));
                if widthoverprime < widthoverdate
                    west(id) = bounds(1);
                    east(id) = bounds(2);
                else
                    west(id) = bounds(2);
                    east(id) = bounds(1);
                end
            elseif crossesprime
                west(id) = min(lon{id}(iseast));
                east(id) = max(lon{id}(iswest));
%                 west(id) = bounds(1);
%                 east(id) = bounds(2);
            elseif crossesdate
                west(id) = bounds(2);
                east(id) = bounds(1);
            end    
            
        end
     
    end
    
    % Figure out west-east bounds of plot that minimize lon span while
    % including all objects
    
    
    for id = 1:ndata
        westbound(id) = west(id);
        eastbounds{id} = [east; west];
        isless = eastbounds{id} < westbound(id);
        eastbounds{id}(isless) = eastbounds{id}(isless) + 360;
        lonspan(id) = max(eastbounds{id}) - westbound(id);
    end
            
    [minval, imin] = min(lonspan);
    
    % Convert bounds to either -180-180 or 0-360 coords as appropriate
    
    westmap = mod(westbound(imin) - In.buffer, 360);
    eastmap = mod(max(eastbounds{imin}) + In.buffer, 360);
    
    if westmap > eastmap
        westmap = westmap - 360;
    end
    
    % Calculate latitude limits
    
    [latmin, latmax] = cellfun(@minmax, lat);
    northmap = max(latmax) + In.buffer;
    southmap = min(latmin) - In.buffer;
             
else
    westmap = In.lonlim(1);
    eastmap = In.lonlim(2);
    southmap = In.latlim(1);
    northmap = In.latlim(2);
end

%-------------------------
% Calculate lon origin and 
% frame limits
%-------------------------
            
if eastmap > westmap
    originlon = (eastmap + westmap)/2;
else
    eastmap = eastmap + 360;
    originlon = mod((eastmap + westmap)/2, 360);
end
    
% Shift limits to make them relative to origin.

flonlim = [westmap eastmap] - originlon;

% Determine frame width.

width = mod(diff(flonlim),360);
width(width == 0) = 360;

% Force flonlim(2) to exceed zero.

west = mod(flonlim(2),360);
west(west == 0) = 360;
flonlim = west + [-width 0];

% But if subtracting 360 would make the limits include or move
% closer to zero, then do so.

if (flonlim(1) > 0) && (360 - flonlim(2) < flonlim(1))
    flonlim = flonlim - 360;
end    

%-------------------------
% Load and trim coastlines
%-------------------------

if strcmp(In.coast, 'm')
    Coast = load('coast');
    coastlat = Coast.lat;
    coastlon = Coast.long;
    [coastlat, coastlon] = maptrimp(coastlat, coastlon, [southmap northmap], [westmap eastmap]);
    
else
    lonlimtemp = [westmap eastmap];
    % lonlimtemp(lonlimtemp > 195) = lonlimtemp(lonlimtemp > 195) - 360;
    lonlimtemp(lonlimtemp > 180) = lonlimtemp(lonlimtemp > 180) - 360;
    if isequal(0, lonlimtemp(1), lonlimtemp(2))
        lonlimtemp = [-180 180];
    end

    if lonlimtemp(2) < lonlimtemp(1)
        Coast1 = gshhs(gshhsfile(In.coast), [southmap northmap], [lonlimtemp(1) 180]);
        Coast2 = gshhs(gshhsfile(In.coast), [southmap northmap], [-180 lonlimtemp(2)]);
        Coast = [Coast1 Coast2];
    else
        Coast = gshhs(gshhsfile(In.coast), [southmap northmap], lonlimtemp);
    end

%     coastlat = {Coast.Lat};
%     coastlon = {Coast.Lon};
%     [coastlat, coastlon] = polyjoin(coastlat, coastlon);
%     [coastlat, coastlon] = maptrimp(coastlat, coastlon, [southmap northmap], [westmap eastmap]);
end

% if In.fillcoast
%     if strcmp(In.coast, 'm')
%         geoshow(coastlat, coastlon, 'displaytype', 'polygon', 'facecolor', 'g');
%     else
%         for ic = 1:length(Coast)
%             geoshow(Coast(ic), 'facecolor', 'g');
%         end
%     end
% %     [f,v] = poly2fv(coastlon, coastlat);
% %     vlat  = v(:,2);
% %     vlon  = v(:,1);
% %     coastlat = vlat(f');
% %     coastlon = vlon(f');
% end    

%-------------------------
% Create map axis
%-------------------------

if isnan(In.axis)
    h.fig = figure;
    axtmp = axes('position', [0 0 1 1]);
else
    axtmp = In.axis;
    axes(axtmp);
end

h.axis = axesm(In.projection, ...
      'MapLatLimit', [southmap northmap], ...
      'MapLonLimit', [westmap eastmap], ...
      'origin', originlon, ...
      'flatlimit', [southmap northmap], ...
      'flonlimit', flonlim, ...
      'parallellabel', 'on', ...
      'meridianlabel', 'on', ...
      'grid', 'on', ...
      'mlinelocation', In.gridint(1), ...
      'plinelocation', In.gridint(2), ...
      'frame', 'on', ...
      In.pv{:});
  
%-------------------------
% Plot data
%-------------------------
  
for il = 1:size(In.line,1)
    h.line(il) = plotm(In.line{il,1}, In.line{il,2});
end

for ic = 1:size(In.contour,1)
    [blah, h.contour(ic)] = contourm(In.contour{ic,1}, In.contour{ic,2}, In.contour{ic,3});
end

for ip = 1:size(In.pcolor,1)
    lattemp = In.pcolor{ip,1};
    lontemp = In.pcolor{ip,2};
    pcotemp = In.pcolor{ip,3};
    
    % Grid if necessary
    
    if ~isequal(size(lattemp), size(lontemp), size(pcotemp))
        [lontemp, lattemp] = meshgrid(lontemp, lattemp);
    end
    
    % Convert to double if necessary
    
    pcotemp = double(pcotemp);
    
    % Squeeze if multiple dimensions
    
    if ndims(pcotemp > 2)
        pcotemp = squeeze(pcotemp);
    end
    
    h.pcolor(ip) = pcolorm(lattemp, lontemp, pcotemp);
end

h.patch = cell(size(In.patch,1),1);
for ip = 1:size(In.patch,1)
    h.patch{ip} = patchesm(In.patch{ip,1}, In.patch{ip,2}, In.patch{ip,3});
    if isfv(ip)
        set(h.patch{ip}, 'edgecolor', 'none');
    end
end

% axis tight;

%-------------------------
% Add coastlines
%-------------------------

landcol = [0.8 0.8 0.8];
if In.fillcoast
    if strcmp(In.coast, 'm')
        h.coast = geoshow(coastlat, coastlon, 'displaytype', 'polygon', 'facecolor', landcol);
    else
        for ic = 1:length(Coast)
            h.coast(ic) = geoshow(Coast(ic), 'facecolor', landcol);
        end
    end
%     h.coast = patchm(coastlat, coastlon, 'facecolor', [.8 .8 .8], 'edgecolor', 'none');
else
    if strcmp(In.coast, 'm')
        h.coast = plotm(coastlat, coastlon, 'color', [.5 .5 .5]);
    else
        for ic = 1:length(Coast)
            h.coast(ic) = plotm(Coast(ic).Lat, Coast(ic).Lon, 'color', [.5 .5 .5]);
        end
    end        
%         
end

%-------------------------
% Subfunctions
%-------------------------

function lon = changecoords(lon)
lon(lon > 180) = lon(lon > 180) - 360;

function b = cells2vec(a)
b = cellfun(@(x) x(:), a, 'uni', 0);
b = catwithnan(b, 1);
