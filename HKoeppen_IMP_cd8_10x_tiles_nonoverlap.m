
function HKoeppen_IMP_cd8_10x_tiles_nonoverlap(varargin)
%function to batch process API calls and image analysis for gSlide scans.
% This script is working on 10x tiles from POPLAR, including creating
% tile-specific outputs in a CSV file.  Overlap is set to 51%, resulting in
% only one tile classification per coordinate, and tile output includes raw
% readout for both CK positive and negative areas

%assigns working directory WDIR according to VARARGIN
for k=1:nargin
    switch varargin{k}
         case '-inputdir'
            if(nargin>=(k+1)) %makes sure there is an input after "-inputdir"
                wdir = varargin{k+1}; %assigns input after "-inputdir" to WDIR
            end
        case '-uid'
            if(nargin>=(k+1)) %makes sure there is an input after "-uid"
                UID = varargin{k+1}; %assigns input after "-uid" to UID
            end
        case '-tilesize'
            if(nargin>=(k+1)) %makes sure there is an input after "-uid"
                tile_dim = varargin{k+1}; %assigns input after "-uid" to UID
            end
        case '-bordersize'
            if(nargin>=(k+1)) %makes sure there is an input after "-uid"
                border_size = varargin{k+1}; %assigns input after "-uid" to UID
            end
        case '-overlapcut'
            if(nargin>=(k+1)) %makes sure there is an input after "-uid"
                overlap_cutoff = varargin{k+1}; %assigns input after "-uid" to UID
            end
        case '-flibyer'
        otherwise
    end
end

slidepath = wdir;
%   textinput = strtrim(textinput);  slidepath = char(textinput(487));  UID = 50;
%   tile_dim = '100' ; border_size = '20';  overlap_cutoff = '.25' ;

%%  Analysis parameters
lowmagpower = 2.5;  %tissue find pass
himagpower = 10;    %CCaspase analysis pass

savefile = strcat(tile_dim , '_' , border_size , '_results_', UID , '.xls'); 
resdir = strcat('/gnet/is2/p01/shares/pathology/microscopy/CALM/CALM_user_data/jea/proj/nanozoomer_matlab/HK_IMP_alltiles_redo_', ...
    tile_dim , '_' , border_size);
heatname = strcat('CD8_CK_map_' , tile_dim, '_' , border_size);
heatname_neg = strcat('CD8_CKNEG_map_' , tile_dim, '_' , border_size);
cd8_ck_maskname = 'CD8_CK_mask';

tile_dim = str2num(tile_dim);
border_size = str2num(border_size); 

% if not using highest mag pixel data, ratio of downsample.  e.g. 20x scan,
% analyzed at 10x = 0.5   Used to calculate overlap ratio
pixel_to_tile_ratio = .5 ; 

% scalar, 0-1, ratio of CK positive/negative area that needs to be present
% to be classified as pos/neg tile
ck_overlap_cutoff = str2num(overlap_cutoff) ;   

mkdir(resdir);
cd(resdir);

mkdir(resdir);
cd(resdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% slide params
splitsep = '<:\|:>';

%% my params

%initialize running totals of areas
tis_area_tot = double(0);
ck_area_tot = double(0);
ckneg_area_tot = double(0);
cd8_ck_area_tot = double(0);
cd8_ckneg_area_tot = double(0);
%dab_int_level = [255 32 28 24 20 16 12 8 4 2 1 -.0001 -.1];
dab_int_level = [2 .20 .15 .12 .08 .06 .04 .02 .01 .005 -.01];

cd8_ck_bins = zeros(size(dab_int_level, 2) - 1, 1) ;
cd8_ckneg_bins = zeros(size(dab_int_level, 2) - 1, 1) ;

ckpos_tile_log = '';
ckneg_tile_log = '';

%clear overlay in gSlide if it exists
opts = struct();
opts.REMOVE = 1;
opts.LAYERNAME = heatname;
OverlayAPI(slidepath,opts);

%clear overlay in gSlide if it exists
opts = struct();
opts.REMOVE = 1;
opts.LAYERNAME = heatname_neg;
OverlayAPI(slidepath,opts);


% %clear overlay in gSlide if it exists
% opts = struct();
% opts.REMOVE = 1;
% opts.LAYERNAME = cd8_ck_maskname;
% OverlayAPI(slidepath,opts);


d30 = strel('disk', 15);
d14 = strel('disk', 7);
d9 = strel('disk', 5);
d3 = strel('square', 3);
d2 = strel('square', 2);

d5 = logical([ 0 1 1 1 0;
    1 1 1 1 1;
    1 1 1 1 1;
    1 1 1 1 1;
    0 1 1 1 0]);

%% Slide info, parse JSON data format into Matlab struct, ASTRUCT

opts = struct();
opts.INFO = 1;
opts.MAGPOWER = lowmagpower;

[status, astruct] = CaptureAPI(slidepath,opts);

%% slide annotations, store attributes in ATTRIB, values in TVALUE
[~, slideannots] = system( sprintf('%s/projectonlayer/getannots.bash -path "%s"',...
   getenv('GSLIDEROOT') , slidepath));
slideannots = regexp(slideannots, splitsep, 'split');
slideannots = reshape(slideannots, [], 1);

[attrib, tvalue] = strtok(slideannots, ':');
for i = 1:size(tvalue, 1)
    tvalue(i,1) = cellstr(tvalue{i,1}(2:end));
end

%% get low-res image

opts = struct();
opts.LEFT = 0;
opts.RIGHT = (astruct.totXpix);
opts.TOP = 0;
opts.BOTTOM = (astruct.totYpix);
opts.MAGPOWER = lowmagpower;

[status, ss_bigmont] = CaptureAPI(slidepath,opts);

%% process low-res image to find tissue

%intensity threshold
tlevel = (230 /255);
tempmont_mask = ~im2bw(ss_bigmont, tlevel);

tlevel = (40 /255);
tempmont_mask = im2bw(ss_bigmont, tlevel) & tempmont_mask;

% %RGB threshold
% redl = 0 ;
% redt = 50 ;
% greenl = 0 ;
% greent = 50 ;
% bluel = 0 ;
% bluet = 50 ;
% del_mask = rgb_thresh(ss_bigmont, redl, redt, greenl, greent, bluel, bluet);
% 
% ss_mask3 = tempmont_mask & ~del_mask ;

ss_mask3 = tempmont_mask;

% density filter to remove non-tumor areas
ss_mask3 = imclose(ss_mask3, d3);
ss_mask3 = imopen(ss_mask3, d3);
ss_mask3 = imclose(ss_mask3, d5);

ss_mask3 = ~(bwareaopen(~ss_mask3, 1000)); % closes small holes in tissue
ss_mask3 = imopen(imclose(ss_mask3, d14), d5);
ss_mask3 = bwareaopen(ss_mask3, 1000); % removes small pieces of tissue
ss_mask3 = imopen(ss_mask3, d5);
ss_mask3 = imclose(ss_mask3, d30);
ss_mask3 = imopen(ss_mask3, d9);
ss_mask3 = bwareaopen(ss_mask3, 1000); % removes small pieces of tissue
ss_mask3 = ~bwareaopen(~ss_mask3, 10000);


%edgemask to bring objects off the edge of the image
emask = true(size(ss_mask3));
emask = bwperim(emask);
ss_mask = ss_mask3 & ~emask;
ss_mask = bwareaopen(ss_mask, 3000); 

ss_mask = bwlabel(im2bw(ss_mask));
stats = regionprops(ss_mask, 'BoundingBox');

if ~isdeployed
    im2 = maskout(ss_bigmont, ss_mask, 3);
end

%%%%%%%%%%%%%
%Manually select SS_MASK region
%  ss_mask3 = roipoly(ss_bigmont);  ss_mask = bwlabel(ss_mask3);  stats = regionprops(ss_mask, 'BoundingBox'); mask_count = 1;

%% create mask of ROIS and keep them


AllRegions = VectorAPI(slidepath, 'GETALL'); %get all ROI ID's to be stamped
if isequal(size(AllRegions, 2) , 1);  AllRegions = {AllRegions} ;  end

del_rois = zeros(size(ss_mask(:,:,1)), 'uint8');  %blank image to be labeled 
add_rois = zeros(size(ss_mask(:,:,1)), 'uint8');  %blank image to be labeled 

for j = 1 : (size(AllRegions, 2))        
        bigmont_mask = gSlide_API_get_ROI_mask(AllRegions{j}.REGIONPTS , astruct.scale , size(del_rois, 1) , size(del_rois, 2)) ; 
        if ( ~isempty(regexpi(AllRegions{j}.NAME, 'Necrosis')) || ...
                ~isempty(regexpi(AllRegions{j}.NAME, 'Exclude')))
            del_rois = bigmont_mask | del_rois;
        elseif ( ( ~isempty(regexpi(AllRegions{j}.NAME, 'Tumor')))  ||...
                ( ~isempty(regexpi(AllRegions{j}.NAME, 'feature'))) )
            add_rois = bigmont_mask | add_rois;
        end
end

add_rois = imresize(add_rois, size(ss_mask), 'nearest');
del_rois = imresize(del_rois, size(ss_mask), 'nearest');

ss_mask = ss_mask & add_rois;
ss_mask = ss_mask & ~del_rois ;
ss_mask = bwareaopen(ss_mask, 1000);
ss_mask = bwlabel(im2bw(ss_mask));
stats = regionprops(ss_mask, 'BoundingBox');

%% process indiviudal tissue pieces
%SS_MASK controlls what portions of the slide are pulled back at himag

%API call to get image data
opts = struct();
opts.MAGPOWER = himagpower;
opts.INFO = 1;
[status, sinfo] = CaptureAPI(slidepath, opts);
AS = (sinfo.umperpixelAtMag) ^ 2 ;   % sq. microns of each pixel

%%
for mask_count = 1 : size(stats, 1)
    
    bigmont_mask = ismember(ss_mask, mask_count);  %make mask with ONLY one region
    hstats = regionprops(bigmont_mask, 'BoundingBox'); %get the bounding box of that region    
    
    %calculate real coordinates for image using BoundingBox
    hstats.t = hstats.BoundingBox(1,2) * astruct.scale ;  %low-mag coordinates multiplied by ration between low-mag and high-mag scale
    hstats.l = hstats.BoundingBox(1,1) * astruct.scale  ;
    hstats.r = (hstats.BoundingBox(1,1) + hstats.BoundingBox(1,3)) * astruct.scale ;
    hstats.b = (hstats.BoundingBox(1,2) + hstats.BoundingBox(1,4)) * astruct.scale ;
    
    opts = struct();
    opts.LEFT = hstats.l;
    opts.RIGHT = hstats.r;
    opts.TOP = hstats.t;
    opts.BOTTOM = hstats.b;
    opts.MAGPOWER = himagpower;
    opts.OUTPUTFORMAT = 'png';
    
    [status, bigmont] = CaptureAPI(slidepath,opts);

    opts.INFO = 1;
    
    [status, sinfo] = CaptureAPI(slidepath, opts);
    AS = (sinfo.umperpixelAtMag) ^ 2 ;   % sq. microns of each pixel
    
    %restrict mask image to ROI specific area
    bigmont_mask = imcrop(bigmont_mask, hstats.BoundingBox);
    bigmont_mask = imresize(bigmont_mask, [size(bigmont, 1) size(bigmont, 2)]);
    
    cd(resdir);

%% find CK positive regions

% Convert RGB image to chosen color space
I = rgb2hsv(bigmont);

% Bright Pink - Define thresholds based on histogram settings
channel1Min = 0.724; channel1Max = 0.893;
channel2Min = 0.275; channel2Max = 0.941;
channel3Min = 0.471; channel3Max = 1.000;

% Create mask based on chosen histogram thresholds
ck_mask = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);

% Deep Purple - Define thresholds based on histogram settings
channel1Min = 0.681; channel1Max = 0.742;
channel2Min = 0.483;channel2Max = 0.818;
channel3Min = 0.659;channel3Max = 0.918;

% Create mask based on chosen histogram thresholds
ck_mask = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max) | ck_mask ;

ck_mask = ck_mask & bigmont_mask ;

ck_mask = imopen(imclose(ck_mask, d2), d2);
ck_mask = imopen(imclose(ck_mask, d5), d2);
ck_mask = imopen(imclose(ck_mask, strel('disk', 7)), d3);
ck_mask = imclose(ck_mask , strel('disk', 11));
ck_mask = ~bwareaopen(~ck_mask, 2000);


    if ~isdeployed        
        im2 = imoverlay(bigmont, bwperim(ck_mask) , [0 1 0]);
        %im2 = imoverlay(im2, bwperim(nmask) , [1 0 0]);
        imshow(im2);
    end
    
    
    %% find CK negative regions
    ckneg_mask = bigmont_mask & ~ck_mask;
    
    
    %% find CD8 positive cells

% Define thresholds for channel 1 based on histogram settings
channel1Min = 0.885;channel1Max = 0.049;
channel2Min = 0.239;channel2Max = 0.865;
channel3Min = 0.165;channel3Max = 0.733;

% Create mask based on chosen histogram thresholds
cd8_mask = ( (I(:,:,1) >= channel1Min) | (I(:,:,1) <= channel1Max) ) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);


% Define thresholds based on histogram settings
channel1Min = 0.926; channel1Max = 0.082;
channel2Min = 0.171;channel2Max = 1.000;
channel3Min = 0.141;channel3Max = 0.682;

% Create mask based on chosen histogram thresholds
cd8_mask = ( (I(:,:,1) >= channel1Min) | (I(:,:,1) <= channel1Max) ) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max) | cd8_mask ;

clear I;

cd8_mask = cd8_mask & bigmont_mask;
        
    if ~isdeployed        
        im2 = imoverlay(bigmont, bwperim(cd8_mask) , [0 1 0]);
        %im2 = imoverlay(im2, bwperim(nmask) , [1 0 0]);
        imshow(im2);
    end

%% calculate CD8 + CK heatmap

    %for a given tile, what is the sum of the CK positive area
    fun = @(block_struct) bwarea(block_struct.data);
    sum_ck = blockproc(ck_mask,[tile_dim tile_dim],fun, 'BorderSize' , [border_size border_size], 'TrimBorder' , false);
    sum_ckneg = blockproc(ckneg_mask,[tile_dim tile_dim],fun, 'BorderSize' , ...
        [border_size border_size], 'TrimBorder' , false );
    
    %for a given tile, what is the sum of the CK positive CD8 area
    ck_mask_d = imdilate(ck_mask, d3); % dilate ck_mask to get adjacent CD8s
    cd8_ck_mask = ck_mask_d & cd8_mask;
    fun = @(block_struct) bwarea(block_struct.data);
    sum_cd8_ck = blockproc(cd8_ck_mask,[tile_dim tile_dim],fun, 'BorderSize' , [border_size border_size], 'TrimBorder' , false);    
    
    cd8_ckneg_mask = ckneg_mask & cd8_mask;
    fun = @(block_struct) bwarea(block_struct.data);
    sum_cd8_ckneg = blockproc(cd8_ckneg_mask,[tile_dim tile_dim],fun, 'BorderSize' , ...
        [border_size border_size], 'TrimBorder' , false );
        
    
    %%  output ALL tile specific measurements
    
    %loop through tiles
    for xi = 1 : size (sum_ck, 2)
        for yi = 1 : size(sum_ck, 1)
            if (sum_cd8_ck(yi, xi)) > 0 || (sum_cd8_ckneg(yi, xi) > 0)
                ckpos_tile_log = sprintf('%s%.0f,%.0f,%f,%.0f,%.0f,%.0f\n' , ckpos_tile_log ,...
                round(hstats.l + (xi * (tile_dim/pixel_to_tile_ratio))) , ...
                round(hstats.t + (yi * (tile_dim/pixel_to_tile_ratio))) , ...
                sum_ck(yi, xi) , sum_cd8_ck(yi, xi)  , ...
                sum_ckneg(yi, xi) , sum_cd8_ckneg(yi, xi) );
            end
        end
    end
    
    
    %%
    
        %find areas with less than ck_overlap_cutoff% contribution of CK mask and remove
    idx = sum_ck < ((tile_dim + (border_size*2)) * (tile_dim + (border_size*2)) * ck_overlap_cutoff);
    sum_ck(idx) = 0;  

    %find ratio of CD8 mask area to CK mask area
    ratio_cd8_ck = sum_cd8_ck ./ sum_ck;
    ratio_cd8_ck(isnan(ratio_cd8_ck)) = -20;  %remove NaN's from dividing 0 by 0
    ratio_cd8_ck(isinf(ratio_cd8_ck)) = -20;  %remove INFs from dividing by 0
    
    %preallocate blankimage
    temp_image = uint8(false(size(ratio_cd8_ck, 1), size(ratio_cd8_ck, 2), 3));

    for dab_cut = 1 : size(dab_int_level, 2) -1
        
        %if isequal(dab_cut, 1); rm = 0; gm = 0; bm = 20; end
        if isequal(dab_cut, 11); rm = 0; gm = 0; bm = 100; end
        if isequal(dab_cut, 10); rm = 0; gm = 0; bm = 250; end
        if isequal(dab_cut, 9); rm = 0; gm = 140; bm = 250; end
        if isequal(dab_cut, 8); rm = 0; gm = 250; bm = 250; end
        if isequal(dab_cut, 7); rm = 140; gm = 250; bm = 140; end
        if isequal(dab_cut, 6); rm = 40; gm = 250; bm = 0; end
        if isequal(dab_cut, 5); rm = 250; gm = 250; bm = 0; end
        if isequal(dab_cut, 4); rm = 250; gm = 140; bm = 0; end
        if isequal(dab_cut, 3);  rm = 250; gm = 0; bm = 0; end
        if isequal(dab_cut, 2);  rm = 250; gm = 0; bm = 250; end
        if isequal(dab_cut, 1);  rm = 250; gm = 250; bm = 250; end
        
        timage = false(size(ratio_cd8_ck));
        timage(ratio_cd8_ck < dab_int_level(dab_cut) & ratio_cd8_ck >= dab_int_level(dab_cut+1) ) = true;
        temp_image = imadd( temp_image , cat(3, uint8(timage).* rm , uint8(timage) .* gm, uint8(timage).* bm));
        
    end
     
    temp_image = imresize(temp_image, tile_dim, 'nearest'); % this will have extra pixels on right and bottom
    temp_image = temp_image(1: size(bigmont, 1) , 1 : size(bigmont, 2) , :) ; % trim extra pixels

    
%%  STATS for binned CK positive percents

   for bin_count = 1 : size(dab_int_level , 2) -1
       cd8_ck_bins(bin_count) = sum(sum( ( ( ratio_cd8_ck < dab_int_level(bin_count)) & ...
           ( ratio_cd8_ck >= dab_int_level(bin_count+1)) ))) + cd8_ck_bins(bin_count);                     
   end
      

        %% add CK POSITIVE  HEATMAP to overlay image in gSlide        
        
        %temp_image = cat(3, uint8(bigmont_mask).* 30, uint8(dab_mask) .*200, uint8(hem_mask).*50);
        temp_image = imresize(temp_image, 0.25 , 'nearest');
        alpha = uint8(rgb2gray(temp_image));
        alpha(alpha>0) = 200;        
        
        opts = struct();
        opts.LAYERNAME = heatname;
        opts.PROJECT = 1;
        opts.SRCIMAGE = temp_image;
        opts.SRCMAG = sinfo.magpower * 0.25;
        %opts.SRCMAG = himagpower ;
        opts.LEFT = hstats.l;
        opts.RIGHT = hstats.r;
        opts.TOP = hstats.t;
        opts.BOTTOM = hstats.b;
        opts.ALPHA = alpha;
        
        [stat, res] = OverlayAPI(slidepath,opts);     
        
        %% calculate CD8 + CK NEGATIVE heatmap


    %for a given tile, what is the sum of the CK negative area
    fun = @(block_struct) bwarea(block_struct.data);
    sum_ckneg = blockproc(ckneg_mask,[tile_dim tile_dim],fun, 'BorderSize' , ...
        [border_size border_size], 'TrimBorder' , false );
    
    %for a given tile, what is the sum of the CK negative CD8 area
    cd8_ckneg_mask = ckneg_mask & cd8_mask;
    fun = @(block_struct) bwarea(block_struct.data);
    sum_cd8_ckneg = blockproc(cd8_ckneg_mask,[tile_dim tile_dim],fun, 'BorderSize' , ...
        [border_size border_size], 'TrimBorder' , false );
    
    
    %find areas with less than ck_overlap_cutoff% contribution of CK mask and remove
    idx = sum_ckneg < ((tile_dim + (border_size*2)) * (tile_dim + (border_size*2)) * ck_overlap_cutoff);
    sum_ckneg(idx) = 0;  
    
        %%  output tile specific measurements
    
    %loop through tiles
    for xi = 1 : size (sum_ckneg, 2)
        for yi = 1 : size(sum_ckneg, 1)
            if (sum_ckneg(yi, xi)) > 0 && (sum_cd8_ckneg(yi, xi) > 0)
            ckneg_tile_log = sprintf('%s%.0f,%.0f,%f,%.0f\n' , ckneg_tile_log ,...
                round(hstats.l + (xi * (tile_dim*2))) , ...
                round(hstats.t + (yi * (tile_dim*2))) , ...
                sum_ckneg(yi, xi) , sum_cd8_ckneg(yi, xi) ) ;
            end
        end
    end

    %%
    %find ratio of CD8 mask area to CK mask area
    ratio_cd8_ckneg = sum_cd8_ckneg ./ sum_ckneg;
    ratio_cd8_ckneg(isnan(ratio_cd8_ckneg)) = -20;  %remove NaN's from dividing 0 by 0
    ratio_cd8_ckneg(isinf(ratio_cd8_ckneg)) = -20;  %remove INFs from dividing by 0
    
    
     %preallocate blankimage
    temp_image = uint8(false(size(ratio_cd8_ckneg, 1), size(ratio_cd8_ckneg, 2), 3));

    for dab_cut = 1 : size(dab_int_level, 2) -1
        
        %if isequal(dab_cut, 13); rm = 0; gm = 0; bm = 20; end
        if isequal(dab_cut, 11); rm = 0; gm = 0; bm = 100; end
        if isequal(dab_cut, 10); rm = 0; gm = 0; bm = 250; end
        if isequal(dab_cut, 9); rm = 0; gm = 140; bm = 250; end
        if isequal(dab_cut, 8); rm = 0; gm = 250; bm = 250; end
        if isequal(dab_cut, 7); rm = 140; gm = 250; bm = 140; end
        if isequal(dab_cut, 6); rm = 40; gm = 250; bm = 0; end
        if isequal(dab_cut, 5); rm = 250; gm = 250; bm = 0; end
        if isequal(dab_cut, 4); rm = 250; gm = 140; bm = 0; end
        if isequal(dab_cut, 3);  rm = 250; gm = 0; bm = 0; end
        if isequal(dab_cut, 2);  rm = 250; gm = 0; bm = 250; end
        if isequal(dab_cut, 1);  rm = 250; gm = 250; bm = 250; end
        
        timage = false(size(ratio_cd8_ckneg));
        timage(ratio_cd8_ckneg < dab_int_level(dab_cut) & ratio_cd8_ckneg >= dab_int_level(dab_cut+1) ) = true;
        temp_image = imadd( temp_image , cat(3, uint8(timage).* rm , uint8(timage) .* gm, uint8(timage).* bm));
        
    end
     
    temp_image = imresize(temp_image, tile_dim, 'nearest'); % this will have extra pixels on right and bottom
    temp_image = temp_image(1: size(bigmont, 1) , 1 : size(bigmont, 2), :) ; % trim extra pixels

        %% add HEATMAP to overlay image in gSlide        
        
        %temp_image = cat(3, uint8(bigmont_mask).* 30, uint8(dab_mask) .*200, uint8(hem_mask).*50);
        temp_image = imresize(temp_image, 0.25, 'nearest');
        alpha = uint8(rgb2gray(temp_image));
        alpha(alpha>0) = 200;        
        
        opts = struct();
        opts.LAYERNAME = heatname_neg;
        opts.PROJECT = 1;
        opts.SRCIMAGE = temp_image;
        opts.SRCMAG = sinfo.magpower * 0.25;
        %opts.SRCMAG = himagpower ;
        opts.LEFT = hstats.l;
        opts.RIGHT = hstats.r;
        opts.TOP = hstats.t;
        opts.BOTTOM = hstats.b;
        opts.ALPHA = alpha;
        
        [stat, res] = OverlayAPI(slidepath,opts); 

        %% add CD8 and CK masks to overlay image in gSlide        
        
        temp_image = cat(3, uint8(cd8_mask).* 200, uint8(ck_mask_d) .*200, uint8(bigmont_mask).*20);
        temp_image = imresize(temp_image, 0.25);
        alpha = uint8(rgb2gray(temp_image));
        alpha(alpha>0) = 200;        
        
        opts = struct();
        opts.LAYERNAME = cd8_ck_maskname ;
        opts.PROJECT = 1;
        opts.SRCIMAGE = temp_image;
        opts.SRCMAG = sinfo.magpower * 0.25;
        %opts.SRCMAG = sinfo.magpower ;
        opts.LEFT = hstats.l;
        opts.RIGHT = hstats.r;
        opts.TOP = hstats.t;
        opts.BOTTOM = hstats.b;
        opts.ALPHA = alpha;
        
        [stat, res] = OverlayAPI(slidepath,opts);     
        
        %% calculate stats        

    %binned histogram CK NEGATIVE
   for bin_count = 1 : size(dab_int_level , 2) -1
       cd8_ckneg_bins(bin_count) = sum(sum( ( ( ratio_cd8_ckneg < dab_int_level(bin_count)) & ...
           ( ratio_cd8_ckneg >= dab_int_level(bin_count+1)) ))) + cd8_ckneg_bins(bin_count);                     
   end
      
   % total pixel areas
   
   ck_area_tot = ck_area_tot + bwarea(ck_mask);
   cd8_ck_area_tot = cd8_ck_area_tot + bwarea(cd8_ck_mask);
   
   ckneg_area_tot = ckneg_area_tot + bwarea(ckneg_mask);
   cd8_ckneg_area_tot = cd8_ckneg_area_tot + bwarea(cd8_ckneg_mask);
   
   tis_area_tot = tis_area_tot + bwarea(bigmont_mask);
   
   
end

%% save slide stats

header = sprintf('^^^^^^^');

loc = strcmp(attrib, 'barcode text');
logger = sprintf('%s', char(tvalue(loc)));

urlOut = char(java.net.URLEncoder.encode(slidepath,'UTF-8'));
encodedpath = strrep(urlOut, '+', '%20');
logger = sprintf('%s^%s^=hyperlink(RC[-1],"gSlide")^', logger, char(strcat('http://gslideviewer.gene.com:3000/#/multiview?PATH1=', encodedpath)));

tval = char(tvalue(strcmp(attrib, 'Group')));
if isempty(tval)
    logger = sprintf('%s^^^^^', logger); else logger = sprintf('%s%s^^^^^', logger, tval); end

tval = char(tvalue(strcmp(attrib, 'Block')));
if ~isempty(tval)
    logger = sprintf('%s^%s', logger, tval); else logger = sprintf('%s^', logger); end
header = sprintf('%s^Block', header); 

tval = char(tvalue(strcmp(attrib, 'Patient Number')));
if ~isempty(tval)
    logger = sprintf('%s^%s', logger, tval); else logger = sprintf('%s^', logger); end
header = sprintf('%s^Patient Number', header); 

tval = char(tvalue(strcmp(attrib, 'Patient.ID')));
if ~isempty(tval)
    logger = sprintf('%s^%s', logger, tval); else logger = sprintf('%s^', logger); end
header = sprintf('%s^Patient.ID', header); 


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

header = sprintf('%s^Tissue area, sq.um^CK area / tissue area * 100^CK NEG area / tissue area * 100^CD8 + CK area / CK area * 100^CD8 + CKNEG area / CKNEG area * 100^CK pos tiles^CK neg tiles', ...
    header);

logger = sprintf('%s^%f^%f^%f^%f^%f^%f^%f', logger , ...
    tis_area_tot * AS , (ck_area_tot / tis_area_tot) * 100 , (ckneg_area_tot / tis_area_tot) * 100 , ...
        (cd8_ck_area_tot / ck_area_tot) * 100 , (cd8_ckneg_area_tot / ckneg_area_tot) * 100 , ...
        sum(cd8_ck_bins) , sum(cd8_ckneg_bins) );

for ii = size(cd8_ck_bins,1) : -1 : 1
    header = sprintf('%s^%s ckpos bin / tot ckpos tiles', header , num2str(dab_int_level(ii)));
end

for ii = size(cd8_ck_bins,1) : -1 : 1
    logger = sprintf('%s^%f', logger , (cd8_ck_bins(ii) / sum(cd8_ck_bins) ) );
end

for ii = size(cd8_ck_bins,1) : -1 : 1
    header = sprintf('%s^%s ckneg bin / tot ckneg tiles', header , num2str(dab_int_level(ii)));
end

for ii = size(cd8_ck_bins,1) : -1 : 1
    logger = sprintf('%s^%f', logger , (cd8_ckneg_bins(ii) / sum(cd8_ckneg_bins)) );
end

header = sprintf('%s\n', header);
logger = sprintf('%s\n', logger);

cd(resdir)
[cvf] = fopen(savefile,'a+');
fwrite(cvf, header);
fwrite(cvf, logger);
fclose(cvf);


%% log tiles
tval = char(tvalue(strcmp(attrib, 'Block')));
sval = char(tvalue(strcmp(attrib, 'Patient Number')));

cvf = fopen( strcat(tval , '_' , sval , '_ckpos_tiles.csv'), 'w');
fwrite(cvf, sprintf('X,Y,ckposarea,cd8_ckposarea,cknegarea,cd8_cknegarea\n') );
fwrite(cvf, ckpos_tile_log);
fclose(cvf);


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Secondary functions
