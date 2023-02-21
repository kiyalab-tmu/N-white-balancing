% Input Image and illuminants must be XYZ color space.
function out = transformNWB(A, illuminants_xyz, illuminants_coordinates, num_matrix, MA)
A_size = size(A);

% calculate the distance between each image pixel and illuminants_coordinates
k = calculateK(cell2mat(illuminants_coordinates), A_size, num_matrix);
 
% reshape A to 2-d array (W*H*3 to (W*H)*3).
A = reshape(A, [A_size(1)*A_size(2), A_size(3)]);
A   = A.'; % (W*H)*3 -> 3*(W*H)

% apply color correction
k = permute(k, [3 2 1]);
destination = MA * whitepoint('d65').';
source = MA * illuminants_xyz.';
weight = (source * k);
out = (MA \ ((destination ./ weight) .* (MA * A)));

% reverse shape of 'out'
out = reshape(out.', [A_size(1), A_size(2), A_size(3)]); % 3*(W*H) -> W*H*3
end

%/////////////////////////////////////////////////////////////////
% calculate coefficient k for S1, S2, ... , Sn.
%/////////////////////////////////////////////////////////////////
function k = calculateK(cordin, A_size, num_matrix)

bw = zeros(A_size(1), A_size(2), num_matrix);
d  = zeros(A_size(1), A_size(2), num_matrix);

for i = 1:num_matrix
    bw(cordin(i,2),cordin(i,1),i) = 1;
    d(:,:,i) = bwdist(bw(:,:,i));
end


k = (1 ./ d) ./ sum(1 ./ d, 3);
k(isnan(k))=1;

k = reshape(k, [1, A_size(1)*A_size(2), num_matrix]);
end
