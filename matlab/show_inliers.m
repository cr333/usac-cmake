%% display inliers for an image pair
%function show_inliers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% paths and constants
%     clear; 
    clc;
    working_dir = 'E:\rraguram\projects\USAC\data\homog\test2';
    im1_fname = 'im1.jpg';
    im2_fname = 'im2.jpg';
    orig_pts_file = 'orig_pts.txt';
    inliers_file = 'inliers.txt';

    skip = 5;   % show subset of inliers based on skip size
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% read in images 
    im1 = imread(fullfile(working_dir, im1_fname));
    im2 = imread(fullfile(working_dir, im2_fname));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% read in original data points
    fid_o = fopen(fullfile(working_dir, orig_pts_file), 'r');
    num_pts = str2num(fgetl(fid_o));
    m1 = zeros(2, num_pts);
    m2 = zeros(2, num_pts);    
    for i = 1:num_pts
        temp = textscan(fgetl(fid_o), '%s');
        m1(1, i) = str2num(temp{1,1}{1});
        m1(2, i) = str2num(temp{1,1}{2});
        m2(1, i) = str2num(temp{1,1}{3});
        m2(2, i) = str2num(temp{1,1}{4});        
    end
    fclose(fid_o);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% read in inlier data
    inliers = textread(fullfile(working_dir, inliers_file));
    inliers_ind = find(inliers > 0);
    n = randperm(length(inliers_ind));
    inliers_ind = inliers_ind(n);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% display images side by side and plot inliers
    I = [];
    padding = 15;
    [M1,N1,K1]=size(im1);
    [M2,N2,K2]=size(im2);
    N3=N1+N2+padding;
    M3=max(M1,M2);
    oj=N1+padding;
    oi=0;
    I(1:M3, 1:N3, 1) = 255; I(1:M3, 1:N3, 2) = 255; I(1:M3, 1:N3, 3) = 255;
    I(1:M1,1:N1,:) = im1;
    I(oi+(1:M2),oj+(1:N2),:) = im2;

    figure(1);
    title_str = sprintf('%s and %s: %d inliers',im1_fname, im2_fname, length(inliers_ind));
    imshow(uint8(I)); set(1,'name', title_str); hold on; axis image; axis off;
    for m = 1:skip+1:length(inliers_ind)
        plot([m1(1,inliers_ind(m)) m2(1,inliers_ind(m))+oj], [m1(2,inliers_ind(m)) m2(2,inliers_ind(m))+oi],'go-','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','g');
    end
    