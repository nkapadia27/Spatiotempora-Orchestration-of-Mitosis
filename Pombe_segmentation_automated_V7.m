%function seg_bac = bacterial_segmentation_automated

filename=uigetfile('.tif','Pick your BF file');
%filename = uigetfile('*_SME_projection.tif');
prompt_dlg = {'H factor', 'Erode Diameter', 'Dilate Diameter', 'Sensitivity Threshold', 'Max Area Threshold','Eccentricity Threshold','Max Intensity Thresh'};
dlg_title = 'Parameters';
dim_box = [1,50];
definput = {'4','1','3','0.65','50000','0.0','1700'};
user_input_seg = inputdlg(prompt_dlg,dlg_title,dim_box,definput);
%std_thresh = str2num(user_input_std_initial{1});
%std_thresh =1800;
%int_thresh = 

h_factor = str2num(user_input_seg{1}); %default 4
strel_rad = str2num(user_input_seg{2});%default 1
strel_rad_dil = str2num(user_input_seg{3});%default 3
ecc_thresh = str2num(user_input_seg{6});%default 0.5
image_read=imread(filename);
max_size = str2num(user_input_seg{5});%default 400
adapt_sens_thresh = str2num(user_input_seg{4});
int_thresh = str2num(user_input_seg{7});%default 0.7
%% 
rd2= image_read;%{1,1}{1,1};
rd2 = uint16(rd2);
 figure, imshow(rd2,[min(min(rd2)), max(max(rd2))]);
%imcontrast
%% 
% se=strel('rectangle',[100,40]);
% tophatFiltered = imtophat(rd2,se);
% %figure, imshow(tophatFiltered)
%imcontrast



%% 
%thresh_crit=adaptthresh(tophatFiltered,0.8,'Statistic','Gaussian');%ddefault 0.02
thresh=imbinarize(rd2,'adaptive', 'sensitivity', adapt_sens_thresh);%'global');%),thresh_crit);
figure,imshow(thresh);
 cell_Bin_stats_thr = regionprops(thresh, rd2, 'MaxIntensity');
 cell_bin_int_thr = cell2mat({cell_Bin_stats_thr.MaxIntensity});
 cell_int_thresh_thr = find(cell_bin_int_thr < int_thresh);
   %I4_thr = imfuse(rd2, uint16(thresh));
   cc_bin_thr = bwconncomp(thresh);
   num_rev_thr = length(cell_int_thresh_thr);
   for ii = 1:num_rev_thr
       idx_remove_thr = cell_int_thresh_thr(ii);
   thresh(cc_bin_thr.PixelIdxList{idx_remove_thr}) = 0;
   end
% % 
% gauss_init = imgaussfilt(uint8(thresh),0.5);
%  thresh_2 = imbinarize (gauss_init,'global');
% % figure, imshow(thresh_2)
%big_thresh = bwareafilt(thresh, [100,50000]);
big_thresh = thresh;
%  stats_big = regionprops(thresh,'Area');
% bin_objects = bwconncomp(big_thresh);
%  numPixels = cellfun(@numel,bin_objects.PixelIdxList);
%  [biggest,idx] = max(numPixels);
% big_thresh(bin_objects.PixelIdxList{idx}) = 0;
%big_thresh = bwareaopen(big_thresh,800);
% big_thresh = ~big_thresh; 
 big_thresh = bwareafilt(big_thresh, [0,max_size]);
figure,imshow(big_thresh)
%% 

big_thresh = imfill(big_thresh,8,'holes');
%bw2 = bwareaopen(big_thresh,300);
se=strel('disk',strel_rad,8);
 %se = strel('rectangle',[15,2]);
bw2 = imerode(big_thresh, se);
% 
% D_thresh = -bwdist(~bw2,'euclidean');
% %figure, imshow(D,[])
% %D(~bw2) = 0;
% I3_thresh = imhmin(D_thresh,h_factor + 8);
% L_thresh = watershed(I3_thresh);
%L(~bw2) = 0;
% 
% bw3_thresh = bw2;
% bw3_thresh(L_thresh == 0) = 0;
% % figure, imshow(bw3_thresh)
% bw2 = bw3_thresh;
%bw2 = imfill(bw2,8,'holes');
figure, imshow (bw2)
%  se=strel('disk',strel_rad,8);
% bw2 = imerode(bw2, se);

%  bw_gauss = imgaussfilt(uint8(bw2),2.0);
%  bw2 = imbinarize (bw_gauss,'global');
%  figure, imshow(bw2)
%bw2 = bwareafilt(bw2, [200, 50000]);
%bw2 = bwareaopen(bw2, 600);
% bw2=imerode(bw2,se);
%  figure, imshow(bw2)
%% 


% bw2 = mat2gray(bw2);
% bw_gauss = imgaussfilt(bw2,2.5);
% figure, imshow(bw_gauss);
% thresh2=imbinarize(bw_gauss,'adaptive','sensitivity',0.4);
% thresh2 = imfill(thresh2, 8,'holes');%),thresh_crit);
% figure, imshow (thresh2)
% bw2 = thresh2;
%% 

%bw2 = bw;

D = -bwdist(~bw2,'euclidean');
figure, imshow(D,[])
%D(~bw2) = 0;
I3 = imhmin(D,h_factor);
L = watershed(I3);
%L(~bw2) = 0;

bw3 = bw2;
bw3(L == 0) = 0;
figure, imshow(bw3)
%figure,imshow(label2rgb(L,'jet','w'))
%bw3=mat2gray(bw2);
%bw4=imgaussfilt(bw3,12.0);%default 1.7
  %figure,imshow(bw4);
%% 

mask = imextendedmin(D,h_factor,8);
imshowpair(bw2,mask,'blend')
D2 = imimposemin(D,mask);
%I4 = imhmin(D2,0.8);
Ld2 = watershed(D2);
bw4 = bw2;
bw4(Ld2 == 0) = 0;
bw4 = bwareafilt(bw4, [0, max_size]);
figure,imshow(bw4)
bw5=bw4;
% bw_gauss = imgaussfilt(uint8(bw4),1.5);
% bw5 = imbinarize (bw_gauss,'global');
%figure, imshow(bw5)
% %% 
se_dil = strel('disk',strel_rad_dil,8);
bw5 = imdilate (bw5,se_dil);
D2 = -bwdist(~bw5,'euclidean');
mask2 = imextendedmin(D2,h_factor,8);
% % % % 
% % % % imshowpair(bw4,mask2,'blend')
D3 = imimposemin(D2,mask2);
 Ld3= watershed(D3);
   bw6 = bw5;
   bw6(Ld3 == 0) = 0;
  
 circ_stats = regionprops(bw6, 'Eccentricity');
 bw6 = bwpropfilt(bw6, 'Eccentricity',[ecc_thresh, 1.0]);
 bw6 = imfill(bw6,8, 'holes');
 bw6 = bwareafilt(bw6, [600, max_size]); 
 cell_Bin_stats = regionprops(bw6, rd2, 'MaxIntensity');
 cell_bin_int = cell2mat({cell_Bin_stats.MaxIntensity});
 cell_int_thresh = find(cell_bin_int < int_thresh);
   I4 = imfuse(rd2, uint16(bw6));
   cc_bin = bwconncomp(bw6);
   num_rev = length(cell_int_thresh);
   for ii = 1:num_rev
       idx_remove = cell_int_thresh(ii);
   bw6(cc_bin.PixelIdxList{idx_remove}) = 0;
   end
   cc_bin2 = imclearborder(bw6,8);
   cc_bin3 = bwconncomp(cc_bin2);
   
   label_bin = labelmatrix(cc_bin3);
   rgb_I4 = label2rgb(label_bin, 'jet','w','shuffle');
  figure,imshow(cc_bin2)
 % I5 = imfuse(rd2, rgb_I4);
  figure, imshow (rgb_I4)
% figure,imshow(bw5)
%% 
% 
% thresh_crit2=adaptthresh(bw4,0.5,'Statistic','Gaussian');%default 0.02
% thresh2=imbinarize(bw4,thresh_crit2);
%  %figure,imshow(thresh2);
% 
% %% 
% bw5=bwareaopen(thresh2,150);%default 100
% bw6=uint16(bw5);
%  %figure,imshow(bw5);
%  %cc=bwconncomp(bw5)
%  stats = regionprops(bw5,rd2,'Image','PixelList','PixelValues');
%  PixelValues = {stats.PixelValues}.';
%  PixelList = {stats.PixelList}.';
% 
% %% 
% 
% image_replace = rd2;
% bin_image = bw5;
% % for i = 1:length(PixelValues)
% %     std_int = std(double(PixelValues{i}));
% %     mean_int = mean(double(PixelValues{i}));
% %     px_coordinates = PixelList{i};
% %     if std_int < std_thresh %| mean_int < int_thresh
% %         for j = 1:length(px_coordinates)
% %         image_replace(px_coordinates(j,2), px_coordinates(j,1)) = 0;
% %         bin_image(px_coordinates(j,2), px_coordinates(j,1)) = 0;
% %         end
% %     else
% %         continue
% %     end
% % end
% bin_image = uint16(bin_image);
% fused_image = imfuse(bin_image, rd2,'blend');
% %figure,imshow(image_replace,[min(min(image_replace)), max(max(image_replace))])
% %figure, imshow(bin_image,[min(min(bin_image)), max(max(bin_image))]) 
% figure, imshow(rd2,[min(min(rd2)), max(max(rd2))]);
% figure, imshow(fused_image)
% savename=strrep(filename,'SME_projection','binary_auto');
% 
% user_input = inputdlg ('Is the segmentation good?','Yeast Segmentation',[1 50],{'NO'});
%     if user_input{1} == 'Y'| user_input {1} == 'y'
%         close all
%         imwrite(bin_image,savename,'tif') 
%         return
%     else 
%         close all
%         return
%         %default_std_thresh = num2str(std_thresh);
%         %user_input_std = inputdlg('Std threshold','Set a new threshold',[1 50],{default_std_thresh});
%         %std_thresh = str2num(user_input_std{1});
%     end
% 
% %imwrite(bin_image,savename,'tif')
% %end
