close all
% constants
CELL_INDEX = 2;
AREA_INDEX = 3;
COL_INDEX = 4;
ROW_INDEX = 5;
MAX_RAD_INDEX = 13;
NUM_PATIENTS = 5;
% patient ids for each cell, whether they come from tumor areas or non-tumor areas
% x and y coordinates of cells and threshold for eliminating too small
% cells
pat_inds = [1792, 1794, 1798, 1843, 1910];
tumor_inds = ['N', 'T', 'N', 'T', 'N'];
chosen_xs = [2341, 2759, 1034, 1192, 62];
chosen_ys = [737, 1341, 2205, 981, 1963];
dummy_radius = 0;
layer = 0;
thresh = 100;
% chosen cells
for i=1:NUM_PATIENTS
    pat_ind = pat_inds(i);
    chosen_x = chosen_xs(i);
    chosen_y = chosen_ys(i);
    dat = importdata(strcat('Nuclei_', num2str(pat_ind), '.txt'));
    dat = dat.data;
    imname = strcat(num2str(pat_ind), '_', tumor_inds(i), '_layer1.tif');
    im = extractCellPatch(imname, chosen_x, chosen_y, dummy_radius, layer, dat);
    pause;
end
% random cells
for i=1:NUM_PATIENTS
    pat_ind = pat_inds(i);
    imname = strcat(num2str(pat_ind), '_', tumor_inds(i), '_layer1.tif');
    dat = importdata(strcat('Nuclei_', num2str(pat_ind), '.txt'));
    dat = dat.data;
    inds = find(dat(:,AREA_INDEX) >= thresh);
    d_sel = dat(inds,:);
    chosenInd = randi(size(d_sel,1));
    realInd = d_sel(chosenInd, CELL_INDEX);
    row = dat(realInd, ROW_INDEX);
    col = dat(realInd, COL_INDEX);
    rad = ceil(dat(realInd, MAX_RAD_INDEX));
    im2 = extractCellPatch(imname, row, col, rad, layer);
    pause;
end