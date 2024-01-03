function [binarized_img] = img_bin(im)

% binarizaing the image with 0.1 threshold value
im_nb = imbinarize(im2gray(im),0.1);

% filling holes in the vessel region
se = strel('disk',1);
im_nb2 = imclose(im_nb,se);

% removing noise after image binarization
im_nb3 = bwareaopen(im_nb2, 100);

% cleaning and smoothing edges of binarized regions
se2 = strel('line',5,180);
im_nb4 = imdilate(im_nb3,se2);


windowSize=10;  % Decide as per your requirements
kernel=ones(windowSize)/windowSize^2;
result=conv2(single(im_nb4),kernel,'same');
result=result>0.5;
im_nb4(~result)=0; 

im_nb5 = bwmorph(im_nb4,'majority');

img_x_lim = size(im_nb5,2);
img_y_lim = size(im_nb5,1);

im_nb5(:,img_x_lim) =  im_nb5(:,img_x_lim-1);
binarized_img = im_nb5;

end