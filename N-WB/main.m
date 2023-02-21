clear
clc
close all

%% options
imgDir = '../IMG/'; % input dir
outDir = '../results_3x3_PCA/'; % output dir
catMethod  = 'bradford'; % select from 'simple', 'vonkries' or 'bradforad'
colorSpace = 'linear-rgb'; % image colorspace
NAADthreshold = 0.0; % set 0 to disable

% illuminant estimation method
% Your custom function must have arguments for an image 'img' and its mask 'mask' if you want to use such a function.
estinFunc = 'illumpca'; % select from 'illumwhite', 'illumgray', 'illumpca', or your custom function for illuminant estimation

blk = [3, 3]; % number of blocks, where prod(blk) = $N$
num_illuminants = prod(blk);

%% process images in imgDir
dirInfo = dir(imgDir);
for idx = 1:length(dirInfo)
    % read struct data
    imgInfo = dirInfo(idx);
    % skip if a read file is dir or a hidden file
    if imgInfo.bytes <= 16384 || imgInfo.isdir
        continue
    end
    
    I_input = imread(strcat(imgInfo.folder, '/', imgInfo.name));
    I_input = im2double(I_input);
    
    sourceColor = cell(num_illuminants,1);
    destinColor = cell(num_illuminants,1);
    sourceColorCoord = cell(num_illuminants,1);
    
    % excute illuminat estimation
    [estimate_illums, illums_coordinates] = blockwiseEstimationV3(I_input, blk(1), blk(2), [], estinFunc, '', NAADthreshold);
    %[estimate_illums, illums_coordinates] = blockwiseEstimationV2(I_input, blk(1), blk(2), estinFunc);
    
    % create cells for source and destination (ground truth) color
    for i = 1:num_illuminants
        sourceColor{i} = estimate_illums(i,:) * 255;                
        sourceColorCoord{i} = illums_coordinates(i,:);
        destinColor{i} = [1 1 1];
    end
    
    % apply N-WB
    illumInfo = {sourceColor{:}; destinColor{:}; sourceColorCoord{:};}.';
    I_output = chromadaptNWB(I_input, illumInfo, 'ColorSpace', colorSpace, 'Method', catMethod);
    
    % save adjusted image
    I_output = xyz2rgb(I_output, 'ColorSpace', colorSpace, 'OutputType', 'uint16');
    if ~isfolder(outDir)
        mkdir(outDir);
    end
    imwrite(I_output , strcat(outDir, imgInfo.name)); % save sRGB image
end
