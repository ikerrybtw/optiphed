function [ cell_im ] = extractCellPatch(image, x, y, radius, layer)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if strcmp(image(end-3:end),'.svs')
   I = imread(image, 'index', layer);
   cell_im = I(y-radius+1:y+radius, x-radius+1:x+radius,:);
   figure;
   imshow(cell_im);
elseif strcmp(image(end-3:end),'.tif')
    if layer ~= 0
        error('Input is tif file, but you specified layer. Please check again.');
    else
        I = imread(image);
        cell_im = I(y-radius+1:y+radius, x-radius+1:x+radius,:);
        figure;
        imshow(cell_im);
    end
else
    error('Unknown file format.');
end
end

