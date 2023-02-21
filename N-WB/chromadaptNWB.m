function  B = chromadaptNWB(varargin)

inputs = parseInputs(varargin{:});

% Normalize illuminants
% This prevents changing the overall brightness of the scene.

illuminants_xyz = cell(inputs.illumCellLen, 2);

for i = 1:inputs.illumCellLen
    % Col1: source white point -----------------------------------------------
    illuminants_xyz{i,1} =  rgb2xyz(inputs.illuminants{i,1}, ...
                                     'ColorSpace', inputs.ColorSpace);
    
    if strcmp(inputs.Normalize, 'True')
    illuminants_xyz{i,1} = illuminants_xyz{i,1} / illuminants_xyz{i,1}(2);
    end
    
    % Col2: destination white point ------------------------------------------
    illuminants_xyz{i,2} =  rgb2xyz(inputs.illuminants{i,2}, ...
                                     'ColorSpace', inputs.ColorSpace);
             
    if strcmp(inputs.Normalize, 'True')
    illuminants_xyz{i,2} = illuminants_xyz{i,2} / illuminants_xyz{i,2}(2);
    end
end

% coordinates ---------------------------------------------------------

illuminants_coordinates = cell(inputs.illumCellLen, 1);

for i = 1:inputs.illumCellLen
    illuminants_coordinates{i,1} = inputs.illuminants{i,3};
end
% illuminants -----------------------------------------------------
illuminants_xyz_mat = cell2mat(illuminants_xyz(:,1));
MA = MakeCAT(inputs.Method);

% color space conversion ----------------------------------------------

A_XYZ = rgb2xyz(inputs.A,             ...
                'ColorSpace', inputs.ColorSpace);
 
% transform -----------------------------------------------------------

B_XYZ = transformNWB(A_XYZ,                ...
                     illuminants_xyz_mat,     ...
                     illuminants_coordinates, ...
                     inputs.illumCellLen,     ...
                     MA);
B     = B_XYZ;
end

%--------------------------------------------------------------------------
function MA = MakeCAT(Method)
% select CAT
if       strcmp(Method, 'bradford')

    MA = [ 0.8951000  0.2664000 -0.1614000; ...
          -0.7502000  1.7135000  0.0367000; ...
           0.0389000 -0.0685000  1.0296000];

elseif  strcmp(Method, 'vonkries')

    MA = [ 0.4002400  0.7076000 -0.0808100; ...
          -0.2263000  1.1653200  0.0457000; ...
           0.0000000  0.0000000  0.9182200];

elseif  strcmp(Method, 'cat02')
    
    MA = [ 0.7328000  0.4296000 -0.1624000; ...
          -0.7036000  1.6975000  0.0061000; ...
           0.0030000  0.0136000  0.9834000];

else

    MA = [ 1.0000000  0.0000000  0.0000000; ...
           0.0000000  1.0000000  0.0000000; ...
           0.0000000  0.0000000  1.0000000];
end
end

%--------------------------------------------------------------------------
function inputs = parseInputs(varargin)

narginchk(2,8);
matlab.images.internal.errorIfgpuArray(varargin{1:2});
parser = inputParser();
parser.FunctionName = mfilename;

% input image
validateImage = @(x) validateattributes(x, ...
    {'single','double','uint8','uint16'}, ...
    {'real','nonsparse','nonempty'}, ...
    mfilename,'A',1);
parser.addRequired('A', validateImage);

%----------------------------------------------------------

% illuminants (col:1-2) + coordinates (col:3)
validateIlluminants = @(x) validateattributes(x, ...
    {'cell'}, ...
    {'ncols', 3}, ...
    mfilename,'illuminants',2);
parser.addRequired('illuminants', validateIlluminants);

%----------------------------------------------------------

validateStringInput = @(x,name) validateattributes(x, ...
    {'char','string'}, ...
    {'scalartext'}, ...
    mfilename, name);

% nameValue 'ColorSpace': 'srgb', 'adobe-rgb-1998' or 'linear-rgb'
validColorSpaces = {'srgb','adobe-rgb-1998','linear-rgb'};
defaultColorSpace = validColorSpaces{1};
validateColorSpace = @(x) validateStringInput(x,'ColorSpace');
parser.addParameter('ColorSpace', ...
    defaultColorSpace, ...
    validateColorSpace);

% nameValue 'Method': 'bradford', 'vonkries' or 'simple'
validMethods = {'simple','vonkries','bradford'};
defaultMethod = validMethods{1};
validateMethod = @(x) validateStringInput(x,'Method');
parser.addParameter('Method', ...
    defaultMethod, ...
    validateMethod);

% nameValue 'Normalize': 'True', 'False'
validNormalize = {'True','False'};
defaultNormalize = validNormalize{1};
validateNormalize = @(x) validateStringInput(x,'Normalize');
parser.addParameter('Normalize', ...
    defaultNormalize, ...
    validateNormalize);

parser.parse(varargin{:});
inputs = parser.Results;

inputs.illumCellLen = size(inputs.illuminants,1);

% additional validation
% 'A' must be a MxNx3 RGB image
validColorImage = (ndims(inputs.A) == 3) && (size(inputs.A,3) == 3);
if ~validColorImage
    error(message('images:validate:invalidRGBImage','A'));
end

% illuminant cannot be black [0 0 0]
for i = 1:inputs.illumCellLen
    if (isequal(inputs.illuminants{i,1}, [0 0 0]) || ...
        isequal(inputs.illuminants{i,2}, [0 0 0]))
    
        error(message('images:awb:illuminantCannotBeBlack'));
   
    end

    
end

inputs.ColorSpace = validatestring( ...
    inputs.ColorSpace, ...
    validColorSpaces, ...
    mfilename, 'ColorSpace');

inputs.Method = validatestring( ...
    inputs.Method, ...
    validMethods, ...
    mfilename, 'Method');

inputs.Normalize = validatestring( ...
    inputs.Normalize, ...
    validNormalize, ...
    mfilename, 'Normalize');
end

%--------------------------------------------------------------------------
