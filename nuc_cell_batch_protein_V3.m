
[filenames_fluor,path_fluorescence] = uigetfile ('.tif','Pick Fluorescent Files','Multiselect','on');
filenames_fluor = filenames_fluor';
[trk_image_files, path_trk_files] = uigetfile('trk-Labelled.tif','Select tracked cell images', 'Multiselect','on');
trk_image_files = trk_image_files';
[nuc_mask_files, path_nuc_mask_files] = uigetfile('.tif','Select nuclear masked images', 'Multiselect','on');
nuc_mask_files = nuc_mask_files';
BD_file = importdata('trk-birth-death.csv');
Div_file = importdata('trk-division.csv');
num_timepoints = length(nuc_mask_files);

max_num_frames = 96; %cut off for timelapse%50 for CDC25, 70 for Cdc13
time_int = 5; %minutes
pixel_size = 0.130;
div_time_cut = 1;   % 1 if you only want to consider cells born in the video                                                                                                                                                                                                                                                                                                                                                                                                    
%div_time_max_cut = 55; 
%length_ext_cut = 1;
per_pixels = 0.15;
med_filter_sz = 3;
length_ext_cut = 1;
%nuc_cell_time = cell(num_time_points_max,1);
nuc_cell_time = {};
BD_data = BD_file.data;
elem_born = find(BD_data(:,2)>div_time_cut);
BD_data_revised = BD_data(elem_born,:);
%elem_split = find(BD_data_revised(:,3)<div_time_max_cut);
%BD_data_revised = BD_data_revised(elem_split,:);
elem_death = find(BD_data_revised(:,3)<max_num_frames);
BD_data_revised_2 = BD_data_revised(elem_death,:);
non_art = find((BD_data_revised_2(:,3)-BD_data_revised_2(:,2)) >0);
BD_data_revised_2 = BD_data_revised_2(non_art,:);
num_cells = size(BD_data_revised_2,1);
div_time = zeros(num_cells,2);
div_length = zeros(num_cells,3);
for i = 1:num_cells
    div_time(i,2) = time_int*(BD_data_revised_2(i,3) - BD_data_revised_2(i,2));
    frame_grab = BD_data_revised_2(i,3) - 1;
    img_select = imread(trk_image_files{frame_grab});
    cell_select = BD_data_revised_2(i,1);
    img_select_rev = img_select;
    img_select_rev(img_select_rev~=cell_select) = 0;
    img_select_rev(img_select_rev == cell_select) = 1;
%     stats = regionprops(img_select_rev, 'MajorAxisLength');
%     div_length(i,3) = stats.MajorAxisLength*pixel_size;
    ferprop_div = bwferet(img_select_rev,'MaxFeretProperties');
    div_length(i,3) = ferprop_div.MaxDiameter*pixel_size;
%     
    frame_grab_born = BD_data_revised_2(i,2);
    img_select_born = imread(trk_image_files{frame_grab_born});
    img_select_born_rev = img_select_born;
    img_select_born_rev(img_select_born_rev~=cell_select) = 0;
    img_select_born_rev(img_select_born_rev == cell_select) =1;
%     stats_born = regionprops(img_select_born_rev, 'MajorAxisLength');
%     div_length(i,2) = stats_born.MajorAxisLength*pixel_size;%length at born
    ferprop_born = bwferet(img_select_born_rev,'MaxFeretProperties');
    div_length(i,2) = ferprop_born.MaxDiameter*pixel_size;%length at born
    
    div_length(i,1) = cell_select;
    div_time(i,1) = cell_select;
end
% div_length_cut_filt = find(div_length(:,3) > div_length_cut);
% div_time = div_time(div_length_cut_filt,:);
% div_length = div_length(div_length_cut_filt,:);
find_nonerrors = find(div_time(:,2) ~=0);
div_time = div_time(find_nonerrors,:);
div_length = div_length(find_nonerrors,:);
%div_length = div_length(find_nonerrors,:);
div_time(:,2) = div_time(:,2)/60;%hours
% find_div_errs = find(div_time(:,2)>1);
% div_time = div_time(find_div_errs,:);
% div_length = div_length(find_div_errs,:);
length_extension = div_length(:,3)-div_length(:,2);
 length_err = find(length_extension > length_ext_cut);
 div_time = div_time(length_err,:);
 div_length = div_length(length_err,:);
 length_extension = length_extension(length_err,:);

% figure,
% histogram(div_time(:,2));
% mean(div_time(:,2));
% histogram(div_length(:,3));
% mean_div = mean(div_length(:,3));
% std_div = std(div_length(:,3));
% cov_value = std_div/mean_div;
% disp(mean_div);
% disp(mean(div_time(:,2)));
% 
% [b1,Sfit] = polyfit(div_length(:,2), length_extension,1);
% [Yfit, delta_fit] = polyconf(b1, div_length(:,2), Sfit);
% mdl = fitlm(div_length(:,2),length_extension);
% yCalc1 = polyval(b1, div_length(:,2));
% figure, scatter(div_length(:,2), length_extension)
% hold on 
% plot(div_length(:,2), yCalc1);
% hold off
%% 
for ii = 1:max_num_frames
    


nuc_mask = imread(strcat(path_nuc_mask_files,nuc_mask_files{ii}));
cell_mask = imread(strcat(path_trk_files,trk_image_files{ii}));
fluor_image = imread(strcat(path_fluorescence,filenames_fluor{ii}));

nuc_stats = regionprops(logical(nuc_mask),fluor_image,'MeanIntensity','PixelValues','Centroid','PixelList','SubarrayIdx','ConvexImage','EquivDiameter','Image');

cell_stats = regionprops(logical(cell_mask),fluor_image,'PixelValues','PixelList','SubarrayIdx','ConvexImage','Image');
cell_stats_mask = regionprops(logical(cell_mask),cell_mask,'PixelValues','PixelList','SubarrayIdx','ConvexImage','Image','MeanIntensity');
num_cells = length(cell_stats(:,1));
num_nuc = length(nuc_stats(:,1));
nuc_array = [];
elem_del = [];
    for i = 1:num_nuc
    centr = round(nuc_stats(i).Centroid);
    nuc_array = [nuc_array;centr];
    end
    nuc_cell_arr = cell(num_cells,5);
    for j = 1:num_cells
        cell_sub = cell_stats(j).SubarrayIdx;
        cell_conv = cell_stats(j).ConvexImage;
        cell_image = fluor_image(cell_sub{:});
        cell_image_1 = times (cell_image, uint16(cell_conv));
        cell_image_nonzero = nonzeros(cell_image_1);
        cell_pixels = cell_stats(j).PixelList;
        cell_pixels_values = cell_stats(j).PixelValues;
        cell_pixels_rev = cell_pixels;
        cell_pixels_values_rev = cell_pixels_values;
        nuc_elements = ismember(nuc_array, cell_pixels,'rows');
        nuc_elements_2 = find(nuc_elements);
        num_nuc_found = length(nuc_elements_2(:,1));
        if num_nuc_found == 0
            
            continue
        end
        nuc_inten_arr = zeros(num_nuc_found,3);
        for i = 1:num_nuc_found
            nuc_conv_img = nuc_stats(nuc_elements_2(i)).ConvexImage;
            nuc_bin_arr = nuc_stats(nuc_elements_2(i)).SubarrayIdx;
            nuc_px_list = nuc_stats(nuc_elements_2(i)).PixelList;
            nuc_px_vals = nuc_stats(nuc_elements_2(i)).PixelValues;
            nuc_px_mean = nuc_stats(nuc_elements_2(i)).MeanIntensity;
            nuc_diam = nuc_stats(nuc_elements_2(i)).EquivDiameter;
            cyto_elements = ismember(cell_pixels_rev, nuc_px_list,'rows');
            cell_pixels_rev(cyto_elements,:) = 0;
            cell_pixels_values_rev(cyto_elements,:) = 0;

    %         nuc_cut_img = fluor_image(nuc_bin_arr{:});
    %         [rowf, colf] = find(nuc_conv_img);
    %         nuc_inten_values_1 = times(nuc_cut_img, uint16(nuc_conv_img));
    %         nuc_inten_values_nonzero = nonzeros(nuc_inten_values_1);
    %         nuc_inten_mean = mean(nuc_inten_values_nonzero);
    %         nuc_inten_mean2 = mean(nuc_px_vals);
            nuc_inten_arr(i,1) = nuc_px_mean;
            nuc_inten_arr(i,2) = length(nonzeros(cell_pixels_values_rev));
            nuc_inten_arr(i,3) = nuc_diam;
        end
        %nuc_cell_arr{j,1} = nuc_elements_2;
        tot_nuc_pixels = sum(nuc_inten_arr(:,2));
        wt_nuc = nuc_inten_arr(:,2)/tot_nuc_pixels;
        wt_nuc_mean = times(wt_nuc,nuc_inten_arr(:,1));
        nuc_cell_arr{j,1} = cell_stats_mask(j).MeanIntensity;
        nuc_cell_arr{j,2} = nuc_inten_arr;
        nuc_cell_arr{j,3} = mean(nuc_inten_arr(:,3)); %ferediameter nucleus
        nuc_cell_arr{j,4} = mean(nonzeros(cell_pixels_values_rev));%mean intensity excluding nucleus
        nuc_cell_arr{j,5} = sum(wt_nuc_mean);
    end
    if isempty(nuc_cell_arr(j,1)) == 1
        continue
    end
    %nuc_cell_time{ii} = nuc_cell_arr;
    nuc_cell_time = [nuc_cell_time;nuc_cell_arr];
end


%% 
emp_row=cellfun('isempty',nuc_cell_time); 
nuc_cell_time(any(emp_row(:,[1,2,3]),2),:)=[];
unique_cells = unique(cell2mat(nuc_cell_time(:,1)));
[ele_find_uniq, pos]  = ismember(unique_cells,div_length(:,1));
unique_cells = unique_cells(ele_find_uniq,:);
num_unique = length(unique_cells);
uniq_cell_arr = cell(num_unique,1);
for jj = 1:num_unique
    uni_cell = unique_cells(jj);
    ele_uni_cell = find (cell2mat(nuc_cell_time(:,1)) == uni_cell);
    uniq_cell_arr{jj} = nuc_cell_time(ele_uni_cell,:);
end
%% 
num_cells_fin = length(uniq_cell_arr(:,1));
 figure
for kk = 1:num_cells_fin
    cell_dat = uniq_cell_arr{kk};
    inten_tot_cell = cell2mat(cell_dat(:,1));
    inten_wt_cell = cell2mat(cell_dat(:,4));
    inten_wt_nuc = cell2mat(cell_dat(:,5));
    time_points =((1:length(inten_wt_nuc))-1);
    
   hold on
    %plot (time_points/max(time_points), inten_wt_cell,time_points/max(time_points), inten_tot_cell, time_points/max(time_points), inten_wt_nuc)
    plot (time_points/max(time_points), inten_wt_cell, time_points/max(time_points), inten_wt_nuc)
end
hold off
%% 
num_cells = length(div_length(:,1));
Fluor_cell = cell(num_cells, 3);%cell_ID, cell_length, cell intensity, nuclear intensity, rate concentration increase, max rate growth


for i = 1:num_cells
    cell_ID = div_length(i,1);
    cell_find = find(BD_data(:,1) == cell_ID);
    cell_born_frame = BD_data(cell_find,2);
    cell_death_frame = BD_data(cell_find,3);
    cell_track = cell_death_frame - cell_born_frame + 1;
    cell_length = [];%zeros(cell_track,1);
    cell_fluor = [];%zeros(cell_track,1);
    nuc_fluor_values = [];%zeros(cell_track,1);
    cell_tot_values = [];
    for j = cell_born_frame:cell_death_frame
      %frame_grab_nuc = imread(strcat(path_nuc_fluor,filenames_nuc_fluor{j}));
      frame_grab_fluor = imread(strcat(path_fluorescence,filenames_fluor{j}));
      %frame_grab_bin = imread(strcat(path_bin,filenames_bin{j}));
      frame_grab_trk_bin = imread(trk_image_files{j});
%       nuc_stats = regionprops(logical(frame_grab_nuc),'Centroid');
%       nuc_centroids = cell2mat({nuc_stats.Centroid}');
      cell_stats = regionprops(logical(frame_grab_trk_bin),frame_grab_trk_bin,'PixelValues', 'PixelList','MeanIntensity','MaxFeretProperties');
      cell_stats_fluor = regionprops(logical(frame_grab_trk_bin),frame_grab_fluor, 'PixelValues','PixelList','MeanIntensity','MaxFeretProperties');
      
      cell_isolate = find(cell2mat({cell_stats.MeanIntensity}')==cell_ID);
      if isempty(cell_isolate)==1
         Fluor_cell{i,1} = [];
          break
      end
      mean_cell_inten = cell2mat({cell_stats_fluor.MeanIntensity}');
      mean_tot_cell_inten = mean_cell_inten(cell_isolate);
      cell_pixels = cell2mat({cell_stats.PixelList}');
      %nuc_find = find((nuc_centroids(:,1) == cell_pixels(cell_isolate,1)) & (nuc_centroids(:,2) == cell_pixels(cell_isolate,2)));
      
      length_vals = cell2mat({cell_stats.MaxFeretDiameter}');%cell length at time 
      
      cell_length = [cell_length;length_vals(cell_isolate)];
      cell_tot_values = [cell_tot_values;mean_tot_cell_inten];

    
      
    end
    if length(cell_length) < (cell_death_frame - cell_born_frame)+1
        continue
    end
    Fluor_cell{i,1} = cell_ID;
    Fluor_cell{i,2} = cell_length*pixel_size;
    Fluor_cell{i,3} = cell_tot_values;
end
Fluor_cell = Fluor_cell(~cellfun(@isempty, Fluor_cell(:,1)), :);
elem_sel = cell2mat(Fluor_cell(:,1));
elem_find_uniq = ismember(unique_cells,elem_sel);
uniq_cell_arr = uniq_cell_arr(elem_find_uniq);
unique_cells = unique_cells(elem_find_uniq );
elem_length = ismember (div_length(:,1), unique_cells);
div_length = div_length(elem_length,:);
div_time = div_time(elem_length,:);
num_cells = length(Fluor_cell(:,1));
figure,
histogram(div_time(:,2));
 mean(div_time(:,2));
 histogram(div_length(:,3));
mean_div = mean(div_length(:,3));
 std_div = std(div_length(:,3));
 cov_value = std_div/mean_div;
disp(mean_div);
 disp(mean(div_time(:,2)))
%% 
diff_array_conc = zeros(num_cells,4);
for ii = 1:num_cells
cell_num = ii;
%med_filter_sz = 2;

cell_dat = uniq_cell_arr{cell_num};
    cell_length_time = (cell2mat(Fluor_cell(cell_num ,2)));
    %cell_length_time = (cell2mat(Fluor_cell(cell_num ,2)));
    nuc_diam_time = cell2mat(cell_dat(:,3))*pixel_size;

    cell_fluor_time = cell2mat(cell_dat(:,4));
    cell_bg_intensitie_mean = min(cell_fluor_time);
    
    nuc_fluor_time = cell2mat(cell_dat(:,5));
    nuc_bg_intensities_mean = min(nuc_fluor_time);
    time_points =((1:length(nuc_fluor_time))-1);
   
    time_vals = time_points*time_int;
    [max_elem,Imax] = max(cell_fluor_time);
    Imax_min = Imax - 1;
    nuc_diam_max = nuc_diam_time(Imax);
    nuc_diam_min = nuc_diam_time (Imax_min);
    nuc_fluor_max = nuc_fluor_time(Imax);
    nuc_fluor_min = nuc_fluor_time(Imax_min);
    cell_fluor_max = cell_fluor_time(Imax);
    cell_fluor_min = cell_fluor_time(Imax_min);
    cell_length_max = cell_length_time(Imax);
    cell_length_min = cell_length_time(Imax_min);
    nuc_abs_min = (4/3)*pi*((nuc_diam_min/2)^3) * (nuc_fluor_min - nuc_bg_intensities_mean);
    nuc_abs_max = (4/3)*pi*((nuc_diam_max/2)^3) * (nuc_fluor_max - nuc_bg_intensities_mean);
    cell_abs_min = ((3.5^2)* cell_length_min*pi - (4/3)*pi*((nuc_diam_min/2)^3))*(cell_fluor_min-cell_bg_intensitie_mean);
    cell_abs_max = ((3.5^2)* cell_length_max*pi - (4/3)*pi*((nuc_diam_max/2)^3))*(cell_fluor_max- cell_bg_intensitie_mean);
    diff_nuc  = nuc_abs_max -nuc_abs_min;
    diff_cyto = cell_abs_max - cell_abs_min;
    conc_nuc_diff = nuc_fluor_max - nuc_fluor_min;
    conc_cyto_diff = cell_fluor_max - cell_fluor_min;
    diff_array_conc(ii,1) = diff_nuc;
    diff_array_conc(ii,2) = diff_cyto;
    diff_array_conc(ii,3) = conc_nuc_diff;
    diff_array_conc(ii,4) = conc_cyto_diff;
    
end
figure, 
scatter(diff_array_conc(:,1), diff_array_conc(:,2))
figure, 
scatter(diff_array_conc(:,3), diff_array_conc(:,4));
%% 

time_int = 5;
cell_num =  23;
pixel_size =0.130;
%med_filter_sz = 2;

cell_dat = uniq_cell_arr{cell_num};
    cell_length_time = (cell2mat(Fluor_cell(cell_num ,2)));
    cell_mean = (cell2mat(Fluor_cell(cell_num ,3)));
    %cell_length_time = (cell2mat(Fluor_cell(cell_num ,2)));
    nuc_diam_time = cell2mat(cell_dat(:,3))*pixel_size;
    nuc_vol_time = (4/3)*pi.*((nuc_diam_time./2).^3);
    cell_vol_time =  ((3.5/2)^2).* cell_length_time.*pi;
    cell_fluor_time = cell2mat(cell_dat(:,4));
    nuc_fluor_time = cell2mat(cell_dat(:,5));
    time_points =((1:length(nuc_fluor_time))-1);
%     cell_fluor_smooth = smoothdata(cell_fluor_time,'movmedian',med_filter_sz);
%     nuc_fluor_smooth = smoothdata(nuc_fluor_time,'movmedian',med_filter_sz);
    time_vals = time_points*time_int;
    [max_elem,Imax] = max(cell_fluor_time(15:end));
    Imax_min = Imax - 1;
    nuc_diam_max = nuc_diam_time(Imax);
    nuc_diam_min = nuc_diam_time (Imax_min);
    nuc_fluor_max = nuc_fluor_time(Imax);
    nuc_fluor_min = nuc_fluor_time(Imax_min);
    cell_fluor_max = cell_fluor_time(Imax);
    cell_fluor_min = cell_fluor_time(Imax_min);
    cell_length_max = cell_length_time(Imax);
    cell_length_min = cell_length_time(Imax_min);
    nuc_bg_intensities_mean = min(nuc_fluor_time);
     cell_bg_intensitie_mean = min(cell_fluor_time);
      nuc_abs_min = (4/3)*pi*((nuc_diam_min/2)^3) * (nuc_fluor_min - nuc_bg_intensities_mean);
    nuc_abs_max = (4/3)*pi*((nuc_diam_max/2)^3) * (nuc_fluor_max - nuc_bg_intensities_mean);
    cell_abs_min = ((3.5^2)* cell_length_min*pi - (4/3)*pi*((nuc_diam_min/2)^3))*(cell_fluor_min-cell_bg_intensitie_mean);
    cell_abs_max = ((3.5^2)* cell_length_max*pi - (4/3)*pi*((nuc_diam_max/2)^3))*(cell_fluor_max- cell_bg_intensitie_mean);
    nuc_abs_min_alt = ((3.5/2)^2)* cell_length_min*pi * 0.08 *nuc_fluor_min;
    nuc_abs_max_alt = ((3.5/2)^2)* cell_length_max*pi * 0.08 *nuc_fluor_max;
    diff_nuc  = nuc_abs_max -nuc_abs_min;
    diff_nuc_alt = nuc_abs_max_alt - nuc_abs_min_alt;
    diff_cyto = cell_abs_max - cell_abs_min;
    conc_nuc_diff = nuc_fluor_max - nuc_fluor_min;
    conc_cyto_diff = cell_fluor_max - cell_fluor_min;
    %[cross_corr,tlags] = crosscorr(nuc_fluor_time, cell_length_time);
    
    
    %figure, stem(tlags*time_int, cross_corr)
  figure,
    hold on
   yyaxis right
   p1 = plot(time_vals/max(time_vals), cell_length_time, 'r-');
    %p2 =   plot(time_vals/max(time_vals), nuc_diam_time, 'ro');
%    hold on 
   %plot(time_vals/max(time_vals),cell_fluor_time./nuc_fluor_time);
   ylabel('(um)')
    yyaxis left
    %p3 = plot(time_vals/max(time_vals), cell_fluor_smooth, 'b-');
    p3 = plot(time_vals/max(time_vals), cell_fluor_time, 'k-');
  p4 = plot(time_vals/max(time_vals), cell_mean,'g-');
  p5 = plot(time_vals/max(time_vals), nuc_fluor_time,'b-');
  %p5 = plot(time_vals/max(time_vals), cell_fluor_time./nuc_fluor_smooth,'g-');
    ylabel('Mean Intensity')
    legend([p1,p3,p4, p5],{'Cell Length','Cytoplasm Intensity', 'Total Mean Cell Intensity','Nucleus Intensity'})
    hold off
%     
   figure, 
   scatter(diff_array_conc(:,1), diff_array_conc(:,2));
   figure, 
   plot(time_vals/max(time_vals), cell_fluor_time./nuc_fluor_time,'b-');
   figure, 
   subplot(3,1,3)
   plot(time_vals/max(time_vals), cell_fluor_time, 'k-');
    title('Cytoplasmic')
    ylabel('Mean Intensity (A.U)')
    xlabel('Time(normalized)')
   hold on 
   subplot(3,1,2)
   plot(time_vals/max(time_vals), nuc_fluor_time, 'k-');
    title('Nuclear')
    ylabel('Mean Intensity (A.U)')
     xlabel('Time(normalized)')
    hold on
    subplot(3,1,1);
    plot(time_vals/max(time_vals), cell_mean, 'k-');
    title('Total Cell')
    ylabel('Mean Intensity (A.U)')
     xlabel('Time(normalized)')