function [ cell_im ] = extractCellPatch(image, row, col, radius, layer, data)
%extractCellPatch: given an image, and the location and size of a cell
%of interest, extracts image patch corresponding to that cell
%   Detailed explanation goes here

% constants
COL_INDEX = 4;
ROW_INDEX = 5;
MAX_RAD_INDEX = 13;
DEFAULT_RADIUS = 160;
CONSTANT_RAD = 5;
LARGE_CONST_RAD = 50;

% for plotting with small or large margin
big_radius = radius + CONSTANT_RAD;
very_big_radius = radius + LARGE_CONST_RAD;
% nargin = 5 means no data matrix, not the best way of coding this but works
if nargin == 5
    if strcmp(image(end-3:end),'.svs')
        I = imread(image, 'index', layer);
        
        cell_im3 = I(row-very_big_radius+1:row+very_big_radius, col-very_big_radius+1:col+very_big_radius, :);
        figure;
        imshow(cell_im3, 'InitialMagnification', round(DEFAULT_RADIUS / very_big_radius) * 100);
        
        cell_im2 = I(row-big_radius+1:row+big_radius, col-big_radius+1:col+big_radius, :);
        figure;
        imshow(cell_im2, 'InitialMagnification', round(DEFAULT_RADIUS / big_radius) * 100);
        
        cell_im = I(row-radius+1:row+radius, col-radius+1:col+radius, :);
        figure;
        imshow(cell_im, 'InitialMagnification', round(DEFAULT_RADIUS / radius) * 100);
    
    elseif strcmp(image(end-3:end),'.tif')
        if layer ~= 0
            error('Input is tif file, but you specified layer. Please check again.');
        else
            I = imread(image);
            
            cell_im3 = I(row-very_big_radius+1:row+very_big_radius, col-very_big_radius+1:col+very_big_radius, :);
            figure;
            imshow(cell_im3, 'InitialMagnification', round(DEFAULT_RADIUS / very_big_radius) * 100);
            
            cell_im2 = I(row-big_radius+1:row+big_radius, col-big_radius+1:col+big_radius, :);
            figure;
            imshow(cell_im2, 'InitialMagnification', round(DEFAULT_RADIUS / big_radius) * 100);
            
            cell_im = I(row-radius+1:row+radius, col-radius+1:col+radius, :);
            figure;
            imshow(cell_im, 'InitialMagnification', round(DEFAULT_RADIUS / radius) * 100);
            
            
        end
    else
        error('Unknown file format.');
    end
end
% this means data was specified
% in this case radius will not be used
if nargin == 6 
    % find closest cell from data matrix
    dists_y = abs(data(:,COL_INDEX) - col);
    dists_x = abs(data(:,ROW_INDEX) - row);
    dist = dists_x + dists_y;
    [~, ind] = min(dist);
    loc_y = data(ind,COL_INDEX);
    loc_x = data(ind,ROW_INDEX);
    rad = data(ind, MAX_RAD_INDEX);
    
    big_radius = rad + CONSTANT_RAD;
    very_big_radius = rad + LARGE_CONST_RAD;
    if strcmp(image(end-3:end),'.svs')
        I = imread(image, 'index', layer);
        
        cell_im3 = I(row-very_big_radius+1:row+very_big_radius, col-very_big_radius+1:col+very_big_radius, :);
        figure;
        imshow(cell_im3, 'InitialMagnification', round(DEFAULT_RADIUS / very_big_radius) * 100);
        
        cell_im2 = I(row-big_radius+1:row+big_radius, col-big_radius+1:col+big_radius, :);
        figure;
        imshow(cell_im2, 'InitialMagnification', round(DEFAULT_RADIUS / big_radius) * 100);
        
        cell_im = I(loc_x-rad+1:loc_x+rad, loc_y-rad+1:loc_y+rad, :);
        figure;
        imshow(cell_im, 'InitialMagnification', round(DEFAULT_RADIUS / rad) * 100);
        
        
    elseif strcmp(image(end-3:end),'.tif')
        if layer ~= 0
            error('Input is tif file, but you specified layer. Please check again.');
        else 
            I = imread(image);
            
            cell_im3 = I(row-very_big_radius+1:row+very_big_radius, col-very_big_radius+1:col+very_big_radius, :);
            figure;
            imshow(cell_im3, 'InitialMagnification', round(DEFAULT_RADIUS / very_big_radius) * 100);
            
            cell_im2 = I(row-big_radius+1:row+big_radius, col-big_radius+1:col+big_radius, :);
            figure;
            imshow(cell_im2, 'InitialMagnification', round(DEFAULT_RADIUS / big_radius) * 100);
            
            cell_im = I(loc_x-rad+1:loc_x+rad, loc_y-rad+1:loc_y+rad, :);
            figure;
            imshow(cell_im, 'InitialMagnification', round(DEFAULT_RADIUS / rad) * 100);
            
        end
    else
        error('Unknown file format.');
    end
end

