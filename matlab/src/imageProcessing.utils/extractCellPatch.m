function [ im2 ] = extractCellPatch(image, row, col, varargin)
%UNTITLED Given image, extracts and displays patch centered around [row, col]
%   Detailed explanation goes here
DEFAULT_RADIUS = 160;
CONSTANT_RAD = 5;
LARGE_CONST_RAD = 50;

default_rad = 10;
default_datatype = '';
default_layer = 1;
default_outline = false;

expected_datatypes = {'CP', 'R', ''};

p = inputParser;
validTxtFile = @(x) strcmp(x(end-3:end), '.txt');
validImage = @(x) strcmp(x(end-3:end), '.svs') || strcmp(x(end-3:end), '.tif');
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
p.addRequired('image', validImage);
p.addRequired('row', validScalarPosNum);
p.addRequired('col', validScalarPosNum);
p.addParameter('data', '.txt',  validTxtFile);
p.addParameter('rad', default_rad, validScalarPosNum);
p.addParameter('datatype', default_datatype, @(x) any(validatestring(x,expected_datatypes)));
p.addParameter('layer', default_layer, validScalarPosNum);
p.addParameter('outline', default_outline, @(x) boolean(x));
parse(p, image, row, col, varargin{:});

image_type = image(end-3:end);
if strcmp(p.Results.datatype, '') && ~strcmp(p.Results.data, '.txt')
    error('Data file given, but data type not specified.')
elseif strcmp(p.Results.datatype, 'CP')
    COL_INDEX = 4;
    ROW_INDEX = 5;
    MAX_RAD_INDEX = 13;
    AREA_INDEX = 3;
    outlines = strcat(image(1:end-4), '_NucleiOutlines.tif');
else
    AREA_INDEX = 2;
    MAX_RAD_INDEX = 7;
    ROW_INDEX = 8;
    COL_INDEX = 9;
    outlines = strcat('Outline_Mask100_HE_',image(1:4), '.png');
end

if strcmp(image_type, '.svs')
    I = imread(image, 'index', p.Results.layer);
else
    I = imread(image);
end

if strcmp(p.Results.data, '.txt') % no data was given
    radius = p.Results.rad;
    big_radius = radius + CONSTANT_RAD;
    very_big_radius = radius + LARGE_CONST_RAD;
    im3 = I(row-very_big_radius+1:row+very_big_radius, col-very_big_radius+1:col+very_big_radius, :);
    im2 = I(row-big_radius+1:row+big_radius, col-big_radius+1:col+big_radius, :);
    im = I(row-radius+1:row+radius, col-radius+1:col+radius, :);
else
    D = importdata(p.Results.data);
    data = D.data;
    % find closest cell, manhattan distance
    dists_y = abs(data(:,COL_INDEX) - col);
    dists_x = abs(data(:,ROW_INDEX) - row);
    dist = dists_x + dists_y;
    [~, ind] = min(dist);
    loc_y = data(ind,COL_INDEX);
    loc_x = data(ind,ROW_INDEX);
    radius = data(ind, MAX_RAD_INDEX);
    big_radius = radius + CONSTANT_RAD;
    very_big_radius = radius + LARGE_CONST_RAD;
    
    im3 = I(row-very_big_radius+1:row+very_big_radius, col-very_big_radius+1:col+very_big_radius, :);
    im2 = I(row-big_radius+1:row+big_radius, col-big_radius+1:col+big_radius, :);
    im = I(loc_x-radius+1:loc_x+radius, loc_y-radius+1:loc_y+radius, :);
end
if p.Results.outline
    outlines = imread(outlines);
    if strcmp(p.Results.datatype, 'R')
        outlines = outlines';
    end
    
    green = zeros(size(im3));
    green(:,:,2) = 1;
    im_outlines = outlines(row-very_big_radius+1:row+very_big_radius, col-very_big_radius+1:col+very_big_radius, :);
    figure;
    imshow(im3, 'InitialMagnification', round(DEFAULT_RADIUS / very_big_radius) * 100);
    hold on;
    h = imshow(green);      set(h, 'AlphaData', im_outlines(:,:,1)*0.8);
    
    green = zeros(size(im2));
    green(:,:,2) = 1;
    im_outlines = outlines(row-big_radius+1:row+big_radius, col-big_radius+1:col+big_radius, :);
    figure;
    imshow(im2, 'InitialMagnification', round(DEFAULT_RADIUS / big_radius) * 100);
    hold on;
    h = imshow(green);      set(h, 'AlphaData', im_outlines(:,:,1)*0.8);
    
    green = zeros(size(im));
    green(:,:,2) = 1;
    im_outlines = outlines(row-radius+1:row+radius, col-radius+1:col+radius, :);
    figure;
    imshow(im, 'InitialMagnification', round(DEFAULT_RADIUS / radius) * 100);
    hold on;
    h = imshow(green);      set(h, 'AlphaData', im_outlines(:,:,1)*0.8);
else
    figure;
    imshow(im3, 'InitialMagnification', round(DEFAULT_RADIUS / very_big_radius) * 100);
    figure;
    imshow(im2, 'InitialMagnification', round(DEFAULT_RADIUS / big_radius) * 100);
    figure;
    imshow(im, 'InitialMagnification', round(DEFAULT_RADIUS / radius) * 100);
end

