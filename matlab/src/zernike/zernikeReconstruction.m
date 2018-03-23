function [ im ] = zernikeReconstruction(data, cellIndex, imname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% reconstruct shape from Zernike
% doesn't work because we only have magnitude info for coefficients
% can be fixed maybe by using measurementtemplate.py as discussed on
% CellProfiler forum
cellData = data(cellIndex, :);
zernikeFeatures = cellData(21:50);

% N = 16;
% N is diameter, necessary for reconstruction of shape
N = 2 * ceil(cellData(13));
im = zeros(N);
% loop over all poly coefficients
for x=0:N-1
    for y=0:N-1
        count = 1;
        % polar coordinates
        r = sqrt((2*x - N + 1)^2 + (2*y - N + 1)^2) / N/2;
        theta = atan((N-1-2*x) / (2*y - N + 1));
        for i=0:9
            for j=0:i
                % odd ones are all 0
                if mod(i-j, 2) == 1
                    continue;
                end
                % r = sqrt((2*x - N + 1)^2 + (2*y - N + 1)^2) / N;
                im(x+1,y+1) = im(x+1,y+1) + zernikeFeatures(count) * zernfun(i, j, r, theta);
                count = count + 1;
            end
        end
    end
end
figure;
imshow(im);
% this is for future, ignore for now
if nargin == 3
    % I = imread(imname);
    % I = rgb2gray(I);
    loc_x = cellData(5);
    loc_y = cellData(4);
    r = ceil(cellData(13));
    im2 = extractCellPatch(imname, loc_x, loc_y, r, 0);
end

end

