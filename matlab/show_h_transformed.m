%% warp and display images using the estimated homography
%function show_h_transformed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% paths and constants
    working_dir = 'E:\rraguram\projects\USAC\data\homog\test1';
    im1_fname = 'im1.jpg';
    im2_fname = 'im2.jpg';
    h_file = 'H.txt'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% read in images 
    im1 = imread(fullfile(working_dir, im1_fname));
    im2 = imread(fullfile(working_dir, im2_fname));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% read in homography
    H = textread(fullfile(working_dir, h_file));
    H = reshape(H(1:9), 3, 3)';
    H = inv(H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% transform image and display

    mytform = maketform('projective', H');
    [im2t xd yd] = imtransform(im2, mytform, 'FillValues', 255);
    figure; imshow(im2t, 'XData', xd, 'YData', yd); hold on;
    h = imshow(im1, gray(256));
    set(h, 'AlphaData', 0.6)
    
    mytform = maketform('projective', inv(H)');
    [im1t xd yd] = imtransform(im1, mytform, 'FillValues', 255);
    figure; imshow(im1t, 'XData', xd, 'YData', yd); hold on;
    h = imshow(im2, gray(256));
    set(h, 'AlphaData', 0.6)    