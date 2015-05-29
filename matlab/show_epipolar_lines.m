%% reads in inliers and f matrix, shows epipolar lines 
%function show_epipolar_lines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% paths and constants
    clear;
    clc;
    close all;
    working_dir = 'E:\rraguram\projects\USAC\data\fundmatrix\test1';
    im1_fname = 'im1.jpg';
    im2_fname = 'im2.jpg';
    orig_pts_file = 'orig_pts.txt';
    inliers_file = 'inliers.txt';
    fund_file = 'F.txt';

    skip = 5;   % show subset of epipolar lines based on skip size

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% read in fundamental matrix
    F = textread(fullfile(working_dir, fund_file));
    F = reshape(F(1:9), 3, 3)';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% show inlier matches and epipolar lines for all images

    % form vectors of inlying data points
    x1 = m1(:,inliers_ind); x1 = [x1; ones(1, length(x1))];
    x2 = m2(:,inliers_ind); x2 = [x2; ones(1, length(x2))];

    % display images side by side
    I = [];
    [M1,N1,K1]=size(im1);
    [M2,N2,K2]=size(im2);
    N3=N1+N2;
    M3=max(M1,M2);
    oj=N1;
    oi=0;
    I(1:M1,1:N1,:) = im1;
    I(oi+(1:M2),oj+(1:N2),:) = im2;

    % step through each matched pair of points and display the
    % corresponding epipolar lines on the two images
    l2 = F*x1;    % epipolar lines in image2
    l1 = F'*x2;   % epipolar lines in image1

   % solve for epipoles
    [U,D,V] = svd(F);
    e1 = hnormalise(V(:,3));
    e2 = hnormalise(U(:,3));

    figure(1); imshow(uint8(im1)); hold on; axis image; axis off;
    for n = 1:skip+1:length(inliers_ind)
         hline(l1(:,n), 'g'); plot(e1(1), e1(2), 'g*');
    end
    figure(2); imshow(uint8(im2)); hold on; axis image; axis off;    
    for n = 1:skip+1:length(inliers_ind)
         hline(l2(:,n), 'g'); plot(e2(1), e2(2), 'g*');
    end   