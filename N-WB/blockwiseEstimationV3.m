function [chrom, pos] = blockwiseEstimationV3(img, mblk, nblk, mask, estinMethod, EstinFuncType, NAADthreshold, CosineSimthreshold, showWarning)

if nargin < 4; mask = []; end
if nargin < 5; estinMethod = 'illumwhite'; end
if nargin < 6; EstinFuncType = 'none'; end
if nargin < 7; NAADthreshold = 0.1; end
if nargin < 8; CosineSimthreshold = 0.9; end
if nargin < 9; showWarning = false; end
% illuminant estimation method
estinIllumsFunc = str2func(estinMethod);

% calc image size
img_size = size(img);

% If the block size is larger than the image size, then throw an error.
assert(img_size(1) >= mblk);
assert(img_size(2) >= nblk);

% number of pixel in one block of input image.
mpix = round(img_size(1) / mblk);
npix = round(img_size(2) / nblk);

% block-wise illuminant estimation
fun = @(block_struct) IllumsWithPosition(block_struct.data, mask, estinIllumsFunc, block_struct.location, EstinFuncType, NAADthreshold, CosineSimthreshold);
out = blockproc(img,[mpix npix],fun);

% remove remaining pixels due to blk_div
if size(out,1) ~= mblk
    out = out(1:mblk, :);
    if showWarning
        disp('some remaining pixels are excluded (row)');
    end
end
if size(out,2) ~= nblk*5
    out = out(:, 1:nblk*5);
    if showWarning
        disp('some remaining pixels are excluded (col)');
    end
end

% for debug
assert(size(out,1) == mblk);
assert(mod(size(out,2), 5) == 0);
assert(size(out,2)/5 == nblk);

% shape of output
chrom = zeros(mblk * nblk, 3);
pos   = zeros(mblk * nblk, 2);
for i_blk = 1:mblk*nblk
    [i, j] = ind2sub([mblk, nblk], i_blk);
    chrom(i_blk,:) = out(i, ((j-1)*5)+1:(j*5)-2);
    pos(i_blk,:)   = [out(i, ((j-1)*5)+4), out(i, ((j-1)*5)+5)];
end
pos(max(isnan(chrom),[], 2),:)=[];
chrom = rmmissing(chrom);
assert(size(chrom,1)==size(pos,1));
end

function out = IllumsWithPosition(img, input_mask, estinIllumsFunc, location, EstinFuncType, NAADthreshold, CosineSimthreshold)
if isempty(input_mask)
    mask = ones(size(img,1), size(img,2));
else
    mask = input_mask(location(1,1):location(1,1)+size(img,1)-1,location(1,2):location(1,2)+size(img,2)-1);
end
NAADs = zeros(size(img,3),1);
for ch = 1:size(img,3)
    NAADs(ch,1) = NAAD(img(:,:,ch));
end
NAADs = NAADs >= NAADthreshold;
NAADs = logical(prod(NAADs));

if NAADs
    img = im2double(img);
    switch EstinFuncType
        case 'returnpos'
            [Illums, pos] = estinIllumsFunc(img,mask);
        otherwise
            switch func2str(estinIllumsFunc)
                case 'illumwhite'
                    Illums = estinIllumsFunc(img, 'Mask', mask);
                case 'illumgray'
                    Illums = estinIllumsFunc(img, 'Mask', mask);
                case 'illumpca'
                    Illums = estinIllumsFunc(img, 'Mask', mask);
                otherwise
                    Illums = estinIllumsFunc(img,mask);
            end
            Illums_1x1x3 = reshape(Illums, [1 1 3]);
    end

    switch EstinFuncType
        case 'centerpos'
            ind_j = round(size(img,1)/2);
            ind_i = round(size(img,2)/2);
        case 'returnpos'
            ind_j = pos(1);
            ind_i = pos(2);
        otherwise
            % calc cosine similarity
            cos_sim = sum(img .* Illums_1x1x3, 3) ./ (vecnorm(img, 2, 3) .* vecnorm(Illums_1x1x3, 2, 3));
            [max_val, max_ind] = max(cos_sim, [], 'all', 'linear');
            if max_val < CosineSimthreshold
                ind_j = round(size(img,1)/2);
                ind_i = round(size(img,2)/2);
            else
                [ind_j, ind_i] = ind2sub(size(cos_sim),max_ind);
            end
    end

    % output
    out = [Illums, (location(1,2)-1)+ind_i, (location(1,1)-1)+ind_j];
else
    out = NaN(1,5);
end
end

function [val, Smax, N] = NAAD(maksed_img_ch)
if size(maksed_img_ch,2) ~= 1
    maksed_img_ch = reshape(maksed_img_ch,[],1);
end
N = numel(maksed_img_ch); % total number of pixel of segment.
S = mean(maksed_img_ch);
Smax = max(maksed_img_ch);
val = (1/(N*S))*sum(abs(maksed_img_ch - S));
end