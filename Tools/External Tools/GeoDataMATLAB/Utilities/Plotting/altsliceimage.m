function handle = altsliceimage(varargin)


%% Input
GDarr = varargin{1};
Coords = varargin{2};
if nargin==3
    alphamaps = varargin{3};
end