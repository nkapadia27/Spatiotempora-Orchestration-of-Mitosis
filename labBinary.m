filenames_imgs = uigetfile ('*.tif', 'Seg Files', 'Multiselect','on');
filenames_imgs = filenames_imgs';
folder_name_labelled = 'Labelled Images';
mkdir(folder_name_labelled);
num_images = length(filenames_imgs);
for i=1:num_images
    img_slice = imread(filenames_imgs{i});
    img_slice_conn = bwconncomp(img_slice, 8);
    lab_img = labelmatrix(img_slice_conn);
    fil_name_lab = sprintf('Labelled_%03d.tif',i);%strrep(filenames_imgs{i},'seg','labelled-0');
    file_save_nm = strcat(folder_name_labelled,'/',fil_name_lab);
    imwrite(lab_img, file_save_nm);
end
 