%% 
[filenames_fluor,path_fluorescence] = uigetfile ('.tif','Pick Fluorescent Files','Multiselect','on');
%[filenames_nuc_fluor,path_nuc_fluor] = uigetfile ('*tif','Pick Nuclear Segmented Files','Multiselect','on');
%[filenames_bin,path_binary] =  uigetfile ('*tif','Pick Binary Files','Multiselect','on');
filenames_fluor = filenames_fluor';
%filenames_nuc_fluor = filenames_nuc_fluor';
%filenames_bin = filenames_bin';
num_images = length(filenames_fluor);
 
BD_file = importdata('trk-birth-death.csv');
%Div_file = importdata('trk-division.csv');
trk_image_files = uigetfile('trk-Labelled.tif','Select tracked images', 'Multiselect','on');
trk_image_files = trk_image_files';
%% 

max_num_frames = 51; %cut off for timelapse%50 for CDC25, 70 for Cdc13
time_int = 10; %minutes
pixel_size = 0.130;
div_time_cut = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ;% 1 if you only want to consider cells born in the video
%div_time_max_cut = 55;
length_ext_cut = 3;
per_pixels = 0.15;
med_filter_sz = 3;

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
    frame_grab = BD_data_revised_2(i,3)-1;%pick cell length 1 frame before (avoids length increases due to force at separation)
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

figure,
histogram(div_time(:,2));
mean(div_time(:,2));
histogram(div_length(:,3));
mean_div = mean(div_length(:,3));
std_div = std(div_length(:,3));
cov_value = std_div/mean_div;
disp(mean_div);
disp(mean(div_time(:,2)));

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

figure
for i = 1:num_cells
    cell_ID = div_length(i,1);
    cell_find = find(BD_data(:,1) == cell_ID);
    cell_born_frame = BD_data(cell_find,2);
    cell_death_frame = BD_data(cell_find,3);
    cell_track = cell_death_frame - cell_born_frame + 1;
    cell_length = [];%zeros(cell_track,1);
    cell_fluor = [];%zeros(cell_track,1);
    cell_inten_mean = [];
    nuc_fluor_values = [];%zeros(cell_track,1);
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
      cell_inten_mean =[cell_inten_mean; cell_intensity_vals(cell_isolate)];
      
    end
    Fluor_cell{i,1} = cell_ID;
    Fluor_cell{i,2} = cell_length*pixel_size;
    Fluor_cell{i,3} = cell_fluor;
    Fluor_cell{i,4} =  nuc_fluor_values;
    Fluor_cell{i,7} = cell_inten_mean;%cell_inten_mean;
    nuc_fluor_smooth = smoothdata(nuc_fluor_values,'movmedian',med_filter_sz);
    cell_fluor_smooth = smoothdata(cell_fluor,'movmedian',med_filter_sz);
    
    time_vals = ((1:length(nuc_fluor_smooth))-1)*time_int;
    
    
    
    hold on
    subplot(1,2,2)
     plot(time_vals/max(time_vals), cell_fluor./nuc_fluor_values,'DisplayName',num2str(cell_length(end)*pixel_size));
    %text(max(time_vals),(nuc_fluor_values./cell_fluor(end)),num2str(cell_length(end)*pixel_size))
    title('Normalized')
    ylabel('C/N Ratio')
    xlabel('Time(normalized)')
    hold on
    subplot(1,2,1)
    plot(time_vals,  cell_fluor./nuc_fluor_values,'DisplayName',num2str(cell_length(end)*pixel_size));
    %text(max(time_vals),(nuc_fluor_values./cell_fluor(end)),num2str(i))
    title('Raw')
    ylabel('C/N Rato')
    xlabel('Time(minutes)')
   


end

hold off
legend show
Fluor_cell = Fluor_cell(~cellfun(@isempty, Fluor_cell(:,1)), :);
num_cells = length(Fluor_cell(:,1));
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
med_filter = 1;
cell_num =3;


    cell_length_time = (cell2mat(Fluor_cell(cell_num ,2)));
 
    cell_fluor_time = cell2mat(Fluor_cell(cell_num ,3));
    nuc_fluor_time = cell2mat(Fluor_cell(cell_num ,4));
    cell_intensities_mean = cell2mat(Fluor_cell(cell_num ,7));
%     cell_fluor_smooth = smoothdata(cell_fluor_time,'movmedian',med_filter_sz);
%     nuc_fluor_smooth = smoothdata(nuc_fluor_time,'movmedian',med_filter_sz);
    time_vals = ((1:length(cell_length_time))-1)*time_int;
    
%     [cross_corr,tlags] = crosscorr(nuc_fluor_time, cell_length_time);
%     
%     
%     figure, stem(tlags*time_int, cross_corr)
%CDK_smooth = smoothdata(cell_fluor_time./nuc_fluor_time,'movmedian',med_filter_sz);
  figure,
    hold on
    yyaxis right
    p1 = plot(time_vals, cell_length_time);
    ylabel('Cell Length (um)')
    yyaxis left
    %p2 = plot(time_vals/max(time_vals), cell_fluor_time./nuc_fluor_time, 'b-');
    p2 = plot(time_vals, cell_fluor_time./nuc_fluor_time, 'b-');
%      p3 = plot(time_vals/max(time_vals), cell_intensities_mean, 'g--');
%     p4 = plot(time_vals, nuc_fluor_smooth,'k-');
%     p5 = plot(time_vals, nuc_fluor_time,'k--');
    ylabel('C/N Ratio')
    legend([p1, p2],{'Cell Length','C/N Ratio'})
    ylim([0.4, 0.8])
    xlim([0, 130])
    
    hold off
    figure, 
    hold on
    yyaxis right
     p6 = plot(time_vals, cell_length_time);
    ylabel('Cell Length (um)')
    yyaxis left
    p7 = plot(time_vals, cell_intensities_mean);
    legend([p6, p7],{'Cell Length','Mean Cell Intensity'})
    hold off
%     %% 
%    save('Fluorescent_data.mat', 'Fluor_cell')
%   save('Div_length.mat', 'div_length');
%    save('Div_time.mat','div_time');

