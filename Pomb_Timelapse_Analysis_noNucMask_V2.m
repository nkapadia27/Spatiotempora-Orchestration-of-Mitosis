%% 
[filenames_fluor,path_fluorescence] = uigetfile ('.tif','Pick Fluorescent Files','Multiselect','on');
%[filenames_nuc_fluor,path_nuc_fluor] = uigetfile ('*tif','Pick Nuclear Segmented Files','Multiselect','on');
%[filenames_bin,path_binary] =  uigetfile ('*tif','Pick Binary Files','Multiselect','on');
filenames_fluor = filenames_fluor';
%filenames_nuc_fluor = filenames_nuc_fluor';
%filenames_bin = filenames_bin';
num_images = length(filenames_fluor);
 
BD_file = importdata('trk-birth-death.csv');
Div_file = importdata('trk-division.csv');
trk_image_files = uigetfile('trk-Labelled.tif','Select tracked images', 'Multiselect','on');
trk_image_files = trk_image_files';
%% 

max_num_frames =96; %50 for CDC25, 70 for Cdc13
time_int = 5; %minutes
pixel_size = 0.130;
div_time_cut = 1;
%div_time_max_cut = 55;
length_ext_cut = 2;
per_pixels = 0.15;


med_filter_sz = 5;

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

[b1,Sfit] = polyfit(div_length(:,2), length_extension,1);
[Yfit, delta_fit] = polyconf(b1, div_length(:,2), Sfit);
mdl = fitlm(div_length(:,2),length_extension);
yCalc1 = polyval(b1, div_length(:,2));
figure, scatter(div_length(:,2), length_extension)
hold on 
plot(div_length(:,2), yCalc1);
hold off

%% 
num_cells = length(div_length(:,1));
Fluor_cell = cell(num_cells, 7);%cell_ID, cell_length, cell intensity, nuclear intensity, rate concentration increase, max rate growth
figure(1)
figure(2)
figure (3)
figure (4)
for i = 1:num_cells
    cell_ID = div_length(i,1);
    cell_find = find(BD_data(:,1) == cell_ID);
    cell_born_frame = BD_data(cell_find,2);
    cell_death_frame = BD_data(cell_find,3);
    cell_track = cell_death_frame - cell_born_frame + 1;
    cell_length = [];%zeros(cell_track,1);
    cell_fluor = [];%zeros(cell_track,1);
    nuc_fluor_values = [];%zeros(cell_track,1);
    cell_tot_mean = [];
    cell_tot_abs = [];
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
      cell_pixels = cell2mat({cell_stats.PixelList}');
      %nuc_find = find((nuc_centroids(:,1) == cell_pixels(cell_isolate,1)) & (nuc_centroids(:,2) == cell_pixels(cell_isolate,2)));
      
      length_vals = cell2mat({cell_stats.MaxFeretDiameter}');%cell length at time 
      cell_intensity_vals = cell2mat({cell_stats_fluor.MeanIntensity}');%mean intensity of cell at time
      
      pix_intensities = {cell_stats_fluor.PixelValues}';
      pix_intensities_cell = pix_intensities{cell_isolate};
     % pix_intensities_iso = cell2mat(pix_intensities_cell);
      pix_sort = sort(pix_intensities_cell(:),'descend');     
      nuc_inten = pix_sort(1:ceil(length(pix_sort )*per_pixels));
      pix_sort_cell = sort(pix_intensities_cell(:),'ascend');     
      cell_inten = pix_sort_cell(1:ceil(length(pix_sort_cell )*(1-per_pixels)));
      
      cell_length = [cell_length;length_vals(cell_isolate)];
      %cell_fluor = [cell_fluor;cell_intensity_vals(cell_isolate)];
      cell_fluor = [cell_fluor;mean(cell_inten)];
      nuc_fluor_values = [nuc_fluor_values;mean(nuc_inten)];
      cell_tot_mean = [cell_tot_mean; mean(pix_intensities_cell)];
      cell_tot_abs = [cell_tot_abs;sum(cell_inten)];
    end
    if length(cell_length) < (cell_death_frame - cell_born_frame)+1
        continue
    end
    Fluor_cell{i,1} = cell_ID;
    Fluor_cell{i,2} = cell_length*pixel_size;
    Fluor_cell{i,3} = cell_fluor;
    Fluor_cell{i,4} =  nuc_fluor_values;
    Fluor_cell{i,5} = cell_tot_abs;
    Fluor_cell{i,7} = cell_tot_mean;
    Cdc13_Cyto_smooth = smoothdata(cell_fluor,'sgolay',3);
    Cdc13_Nuc_smooth = smoothdata(nuc_fluor_values ,'sgolay',3);
   Cdc13_total_cell_smooth = smoothdata( cell_tot_mean,'sgolay',3);
%     rate_conc = diff(nuc_fluor_smooth)/time_int;
%     max_rate_conc = max(rate_conc);
%     rate_max_elem = find(rate_conc==max_rate_conc);
%     
%     rate_max_elem = rate_max_elem (1) + 1;%in case there's multiple, pick the first one
%     Fluor_cell{i,5} = max_rate_conc;
%     rate_growth = diff(cell_length*pixel_size)/time_int;
%     max_rate_growth = max(rate_growth);
%     Fluor_cell{i,6} = max_rate_growth;
   time_vals = ((1:length(Cdc13_Nuc_smooth))-1)*time_int;
     norm_time = time_vals/max(time_vals);
%     time_idx = [time_vals(rate_max_elem), max(time_vals)];
%     plot_idx = [nuc_fluor_smooth(rate_max_elem),nuc_fluor_smooth(end)];
%     res_show = {strcat('Max:', num2str(max_rate_conc)), num2str(i)};
    figure(1)
    hold on
    subplot(2,1,1)
     plot(time_vals, nuc_fluor_values)%, time_vals,nuc_fluor_values);
    %text(max(time_vals),(nuc_fluor_smooth(end)),num2str(i))
    title('Nucleus')
    ylabel('Mean Intensity (A.U)')
    xlabel('Time(minutes)')
    hold on
    subplot(2,1,2)
    plot(time_vals, cell_fluor)%, time_vals, cell_fluor);
    %text(max(time_vals),(cell_fluor_smooth(end)),num2str(i))
    title('Cytoplasm')
    ylabel('Mean Intensity (A.U)')
    xlabel('Time(minutes)')
    figure(2)
    hold on
    subplot(2,1,1)
     plot(norm_time,Cdc13_Nuc_smooth)%, time_vals/max(time_vals),nuc_fluor_values);
    %text(max(time_vals),(nuc_fluor_smooth(end)),num2str(i))
    title('Nucleus')
    ylabel('Mean Intensity (A.U)')
    xlabel('Time(normalized)')
    hold on
    subplot(2,1,2)
    plot(norm_time, Cdc13_Cyto_smooth)%, time_vals/max(time_vals), cell_fluor);
    %text(max(time_vals),(cell_fluor_smooth(end)),num2str(i))
    title('Cytoplasm')
    ylabel('Mean Intensity (A.U)')
    xlabel('Time(normalized)')
    figure(3), 
    plot(norm_time,Cdc13_Nuc_smooth)
    hold on
    figure (4)
    plot(norm_time, Cdc13_Cyto_smooth)
    hold on
end
hold off

 Fluor_cell = Fluor_cell(~cellfun(@isempty, Fluor_cell(:,1)), :);
 elem_sel = cell2mat(Fluor_cell(:,1));
 elem_find_uniq = ismember(div_length(:,1),elem_sel);
 div_length = div_length(elem_find_uniq,:);
 div_time = div_time(elem_find_uniq,:);
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
% figure
% for i=1:num_cells
%     
%     %cell_length_time = (cell2mat(Fluor_cell(i ,2)))*pixel_size;
%  
%     cell_fluor_time = cell2mat(Fluor_cell(i ,4));
%      cell_fluor_smooth = smoothdata(cell_fluor_time,'movmedian',med_filter_sz);
%     time_vals = (1:length(cell_fluor_time))*time_int;
%   
%     hold on
%     %yyaxis right
%     %plot(time_vals, cell_length_time);
%     %ylabel('Cell Length (um)')
%     %yyaxis left
%     %plot(time_vals, cell_fluor_time)%
%     plot(time_vals, cell_fluor_smooth);
%     text(max(time_vals),(cell_fluor_time(end)),num2str(i))
%     ylabel('Cell Mean Intensity (A.U)')
% end
% hold off
%% 

cell_num = 6;
time_int = 5; %minutes
pixel_size = 0.130;
    cell_length_time = (cell2mat(Fluor_cell(cell_num ,2)));
    cell_fluor_time = cell2mat(Fluor_cell(cell_num ,3));
    nuc_fluor_time = cell2mat(Fluor_cell(cell_num ,4));
    cell_tot_time = cell2mat(Fluor_cell(cell_num,7));
    cell_tot_abs_time = cell2mat(Fluor_cell(cell_num,5));
   Cdc13_Cyto_smooth = smoothdata(cell_fluor_time,'sgolay',3);
    Cdc13_Nuc_smooth = smoothdata(nuc_fluor_time ,'sgolay',3);
   Cdc13_total_cell_smooth = smoothdata(cell_tot_time,'sgolay',3);
    time_vals = ((1:length(cell_length_time))-1)*time_int;
    
    %[cross_corr,tlags] = crosscorr(nuc_fluor_time, cell_length_time);
    
    
    %figure, stem(tlags*time_int, cross_corr)
  figure,
    hold on
    yyaxis right
    p1 = plot(time_vals/max(time_vals), cell_length_time,'Color',[0.6350 0.0780 0.1840],'LineStyle','-','LineWidth',1.0);
    ylabel('Cell Length (um)')
    yyaxis left
    %p2 = plot(time_vals, cell_fluor_smooth, 'b-');
    p3 = plot(time_vals/max(time_vals), cell_fluor_time, 'Color',[0 0.4470 0.7410],'LineStyle','-','LineWidth',1.0);
    p4 = plot(time_vals/max(time_vals), cell_tot_time,'Color',[0.4660 0.6740 0.1880],'LineStyle','-','LineWidth',1.0);
    p5 = plot(time_vals/max(time_vals), nuc_fluor_time,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0);
    ylabel('Mean Intensity (A.U)')
    legend([p1, p3, p4, p5],{'Cell Length','Cytoplasm Intensity', 'Mean Cell Intensity', 'Nucleus Intensity'})
    hold off
    figure,
    subplot(3,1,2)
    plot(time_vals/max(time_vals), Cdc13_Nuc_smooth,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0);
    title('Nucleus')
    hold on
    subplot(3,1,1)
    plot(time_vals/max(time_vals), Cdc13_total_cell_smooth,'Color',[0.4660 0.6740 0.1880],'LineStyle','-','LineWidth',1.0);
    title('Total Mean Cell Intensity')
    hold on
%     subplot(4,1,2)
%     plot(time_vals/max(time_vals), cell_tot_abs_time,'Color',[0.6350 0.0780 0.1840],'LineStyle','-','LineWidth',1.0);
%     title('Total Absolute Intensity')
%     hold on
    subplot(3,1,3)
    plot(time_vals/max(time_vals), Cdc13_Cyto_smooth, 'Color',[0 0.4470 0.7410],'LineStyle','-','LineWidth',1.0);
    title('Cytoplasm')
    hold off
    figure, 
    plot(time_vals/max(time_vals), Cdc13_Cyto_smooth)
    figure, 
    plot(time_vals/max(time_vals), Cdc13_Nuc_smooth)
    
figure, plot (Cdc13_total_cell_smooth)
    %% 
%   save('Fluorescent_data.mat', 'Fluor_cell')
%  save('Div_length.mat', 'div_length');
%   save('Div_time.mat','div_time');

