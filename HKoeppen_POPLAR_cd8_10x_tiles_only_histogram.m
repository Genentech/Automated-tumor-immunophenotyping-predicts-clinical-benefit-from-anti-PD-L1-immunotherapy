function HKoeppen_POPLAR_cd8_10x_tiles_only_histogram(varargin)
%function to batch process API calls and image analysis for gSlide scans.
% This script is working on 10x tiles from POPLAR, including creating
% tile-specific outputs in a CSV file.  Writes out tiles for DL train/test

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
        case '-overlapcut'
            if(nargin>=(k+1)) %makes sure there is an input after "-uid"
                overlap_cutoff = varargin{k+1}; %assigns input after "-uid" to UID
            end
        otherwise
    end
end

slidepath = wdir;
%   textinput = strtrim(textinput);  slidepath = char(textinput(4));  UID = 50;
%   tile_dim = '100'; overlap_cutoff = '.25';

%%  Analysis parameters
lowmagpower = 2.5;  %tissue find pass
himagpower = 10;    %CCaspase analysis pass



%directory to save tiled images
resdir = strcat('/gne/data/pathology/CALM/CALM_user_data/jea/proj/nanozoomer_matlab/HK_POPLAR_alltiles_redo_', ...
    tile_dim );

%% slide annotations, store attributes in ATTRIB, values in TVALUE

splitsep = '<:\|:>';
[~, slideannots] = system( sprintf('%s/projectonlayer/getannots.bash -path "%s"',...
   getenv('GSLIDEROOT') , slidepath));
slideannots = regexp(slideannots, splitsep, 'split');
slideannots = reshape(slideannots, [], 1);

[attrib, tvalue] = strtok(slideannots, ':');
for i = 1:size(tvalue, 1)
    tvalue(i,1) = cellstr(tvalue{i,1}(2:end));
end


% POPLAR and OAK slide names
[~, slidename] = fileparts(slidepath) ; 

scratchdir = strcat('/gstore/scratch/u/jea/HKoeppen_POPLAR_MIL_stats/' , slidename ) ;

% file to contain histogram counts
savefile = strcat(slidename , '_', tile_dim , '_histstats_', UID , '.mat'); 

%cmd = sprintf('rm -rf %s', scratchdir);
%system(cmd);

%%
if ~exist(scratchdir, 'dir')
mkdir(scratchdir)


% scalar, 0-1, ratio of CK positive/negative area that needs to be present
% to be classified as pos/neg tile
ck_overlap_cutoff = str2num(overlap_cutoff) ;  
tilesize = str2num(tile_dim);

mkdir(resdir);
cd(resdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% my params

%initialize running totals of areas
r_counts = zeros(1, 256, 'double');
g_counts = zeros(1, 256, 'double');
b_counts = zeros(1, 256, 'double');


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
if isstruct(AllRegions)
    outRegions = AllRegions(1);
   for j = 2: size( AllRegions, 2)
       outRegions = {outRegions , AllRegions(j)};
   end
   AllRegions = outRegions;
end

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
    

    % create and save tiles for deep learning
          bigmont_mask_pad = padarray(bigmont_mask, [tilesize tilesize] , 'post' );
          fun = @(block_struct) process_block_mask_pct(block_struct, bigmont_mask_pad,...
              scratchdir, slidename , mask_count , 'tif' , ck_overlap_cutoff);
    
          blockproc(bigmont, [tilesize tilesize] , fun , ...
              'PadPartial' , true , 'USeParallel' , false , 'DisplayWaitbar' , true, ...
              'PadMethod' , 'symmetric');

end

%%
cmd = sprintf('chmod 777 -R %s', scratchdir);

system(cmd);
%% Create histogram counts

    cd(resdir);

% get list of tile names
tile_list = file_list(scratchdir, 'tif', 0);

% loop through each tile, add counts to histogram

for tidx = 1 : size(tile_list, 2)
    imgin = imread(tile_list{tidx});

    r_counts = r_counts + histcounts(imgin(:,:,1), 256);
    g_counts = g_counts + histcounts(imgin(:,:,2), 256);
    b_counts = b_counts + histcounts(imgin(:,:,3), 256);
end


% save histogram
save(savefile , "b_counts" , "g_counts" , "r_counts");


end
