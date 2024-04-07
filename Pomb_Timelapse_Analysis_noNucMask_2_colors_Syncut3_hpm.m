%% 
[filenames_fluor,path_fluorescence] = uigetfile ('.tif','Pick SynCut Fluorescent Files','Multiselect','on');
[filenames_Prot2_fluor,path_Prot2_fluor] = uigetfile ('.tif','Pick Cdc13 Fluorescent Files','Multiselect','on');
%[filenames_bin,path_binary] =  uigetfile ('*tif','Pick Binary Files','Multiselect','on');
filenames_fluor = filenames_fluor';
filenames_Prot2_fluor = filenames_Prot2_fluor';
%filenames_bin = filenames_bin';
num_images = length(filenames_fluor);
 
BD_file = importdata('trk-birth-death.csv');
Div_file = importdata('trk-division.csv');
trk_image_files = uigetfile('trk-Labelled*.tif','Select tracked images', 'Multiselect','on');
trk_image_files = trk_image_files';
%% 

max_num_frames =75;%120 usually
time_int = 10; %minutes
pixel_size = 0.130;
div_time_cut = 1;%usually 1 if you only want cells born in the video
length_ext_cut = 3;
per_pixels = 0.15;

BD_data = BD_file.data;
elem_born = find(BD_data(:,2)>div_time_cut);
BD_data_revised = BD_data(elem_born,:);
elem_death = find(BD_data_revised(:,3)<max_num_frames);
BD_data_revised_2 = BD_data_revised(elem_death,:);
num_cells = size(BD_data_revised_2,1);
div_time = zeros(num_cells,2);
div_length = zeros(num_cells,3);
for i = 1:num_cells
    div_time(i,2) = time_int*(BD_data_revised_2(i,3) - BD_data_revised_2(i,2));
    frame_grab = BD_data_revised_2(i,3);
    img_select = imread(trk_image_files{frame_grab});
    cell_select = BD_data_revised_2(i,1);
    img_select_rev = img_select;
    img_select_rev(img_select_rev~=cell_select) = 0;
    img_select_rev(img_select_rev == cell_select) = 1;
    %stats = regionprops(img_select_rev, 'MajorAxisLength');
    %div_length(i,3) = stats.MajorAxisLength*pixel_size;
    ferprop_div = bwferet(img_select_rev,'MaxFeretProperties');
    div_length(i,3) = ferprop_div.MaxDiameter*pixel_size;
    
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
length_extension = div_length(:,3)-div_length(:,2);
 length_err = find(length_extension > length_ext_cut);
div_time = div_time(length_err,:);
div_length = div_length(length_err,:);
 length_extension = length_extension(length_err,:);
% length_extension = div_length(:,3)-div_length(:,2);
% length_err = find(length_extension > length_ext_cut);
% div_time = div_time(length_err,:);
% div_length = div_length(length_err,:);
% length_extension = length_extension(length_err,:);

figure (4),

histogram(div_time(:,2));
mean(div_time(:,2));
figure(6)
histogram(div_length(:,3));
mean_div = mean(div_length(:,3));
std_div = std(div_length(:,3));
cov_value = std_div/mean_div;
disp(mean_div);
disp(mean(div_time(:,2)));


[b1,Sfit] = polyfit(div_length(:,2), length_extension,1);
[Yfit, delta_fit] = polyconf(b1, div_length(:,2), Sfit);
yCalc1 = polyval(b1, div_length(:,2));
figure(5), 
scatter(div_length(:,2), length_extension)
hold on 
plot(div_length(:,2), yCalc1);
hold off

%% 
num_cells = length(div_length(:,1));
Fluor_cell = cell(num_cells, 7); %cell_ID, cell_length, cell intensity, nuclear intensity
for i = 1:num_cells
    cell_ID = div_length(i,1);
    cell_find = find(BD_data(:,1) == cell_ID);
    cell_born_frame = BD_data(cell_find,2);
    cell_death_frame = BD_data(cell_find,3);
    if cell_death_frame ~= 80
        continue
    end
    cell_track = cell_death_frame - cell_born_frame + 1;
    cell_length = [];%zeros(cell_track,1);
    cell_fluor = [];%zeros(cell_track,1);
    nuc_fluor_values = [];%zeros(cell_track,1);
    cell_fluor_2 = [];
    nuc_fluor_values_2 = [];
    cell_fluor_2_total = [];
    for j = cell_born_frame:cell_death_frame
      frame_grab_Prot2 = imread(strcat(path_Prot2_fluor,filenames_Prot2_fluor{j}));
      frame_grab_fluor = imread(strcat(path_fluorescence,filenames_fluor{j}));
      %frame_grab_bin = imread(strcat(path_bin,filenames_bin{j}));
      frame_grab_trk_bin = imread(trk_image_files{j});
%       nuc_stats = regionprops(logical(frame_grab_nuc),'Centroid');
%       nuc_centroids = cell2mat({nuc_stats.Centroid}');
      cell_stats = regionprops(logical(frame_grab_trk_bin),frame_grab_trk_bin, 'PixelList','MeanIntensity','MaxFeretPropertie');
      cell_stats_fluor = regionprops(logical(frame_grab_trk_bin),frame_grab_fluor, 'PixelValues','PixelList','MeanIntensity','MaxFeretPropertie');
      cell_stats_fluor_2 = regionprops(logical(frame_grab_trk_bin),frame_grab_Prot2, 'PixelValues','PixelList','MeanIntensity','MaxFeretPropertie');
      
      cell_isolate = find(cell2mat({cell_stats.MeanIntensity}')==cell_ID);
      if isempty(cell_isolate)==1
         Fluor_cell{i,1} = [];
          break
      end
      cell_pixels = cell2mat({cell_stats.PixelList}');
      %nuc_find = find((nuc_centroids(:,1) == cell_pixels(cell_isolate,1)) & (nuc_centroids(:,2) == cell_pixels(cell_isolate,2)));
      
      length_vals = cell2mat({cell_stats.MaxFeretDiameter}');%cell length at time 
      %cell_intensity_vals = cell2mat({cell_stats_fluor.MeanIntensity}');%mean intensity of cell at time
  %Colour 1    
      pix_intensities = {cell_stats_fluor.PixelValues}';
      pix_intensities_cell = pix_intensities{cell_isolate};
     % pix_intensities_iso = cell2mat(pix_intensities_cell);
      pix_sort = sort(pix_intensities_cell(:),'descend');     
      nuc_inten = pix_sort(1:ceil(length(pix_sort )*per_pixels));
      pix_sort_cell = sort(pix_intensities_cell(:),'ascend');     
      cell_inten = pix_sort_cell(1:ceil(length(pix_sort_cell )*(1-per_pixels)));
  %Colour 2
      pix_intensities_2 = {cell_stats_fluor_2.PixelValues}';
      pix_intensities_cell_2 = pix_intensities_2{cell_isolate};
     % pix_intensities_iso = cell2mat(pix_intensities_cell);
      pix_sort_2 = sort(pix_intensities_cell_2(:),'descend');     
      nuc_inten_2 = pix_sort_2(1:ceil(length(pix_sort_2 )*per_pixels));
      pix_sort_cell_2 = sort(pix_intensities_cell_2(:),'ascend');     
      cell_inten_2 = pix_sort_cell_2(1:ceil(length(pix_sort_cell_2 )*(1-per_pixels)));
      
      
      cell_length = [cell_length;length_vals(cell_isolate)];
      %cell_fluor = [cell_fluor;cell_intensity_vals(cell_isolate)];
      cell_fluor = [cell_fluor;mean(cell_inten)];
      nuc_fluor_values = [nuc_fluor_values;mean(nuc_inten)];
       cell_fluor_2 = [cell_fluor_2;mean(cell_inten_2)];
      nuc_fluor_values_2 = [nuc_fluor_values_2;mean(nuc_inten_2)];
      cell_fluor_2_total = [cell_fluor_2_total;mean(pix_intensities_cell_2)];
      
    end
    if length(cell_length) < (cell_death_frame - cell_born_frame)+1
        continue
    end
    Fluor_cell{i,1} = cell_ID;
    Fluor_cell{i,2} = cell_length;
    Fluor_cell{i,3} = cell_fluor;
    Fluor_cell{i,4} =  nuc_fluor_values;
    Fluor_cell{i,5} = cell_fluor_2;
    Fluor_cell{i,6} = nuc_fluor_values_2;
    Fluor_cell{i,7} = cell_fluor_2_total;
    
end
Fluor_cell = Fluor_cell(~cellfun(@isempty, Fluor_cell(:,1)), :);
num_cells = length(Fluor_cell(:,1));
%% Protein 1(Syncut)
rate_arr = zeros(num_cells,3);
figure(1)
figure(2)
figure(3)
for i=1:num_cells
    
    cell_length_time = (cell2mat(Fluor_cell(i ,2)))*pixel_size;
    
    cell_fluor_time = cell2mat(Fluor_cell(i ,3));
    nuc_fluor_time = cell2mat(Fluor_cell(i,4));
    time_vals = ((1:length(cell_fluor_time))-1)*time_int;
    CDK_activity = nuc_fluor_time./cell_fluor_time;
    rate_CDK = diff(CDK_activity)./diff(time_vals)';
    rate_CDK = round(rate_CDK,3);
    max_rate = max(rate_CDK);
    max_elem = find(rate_CDK == max_rate);
    max_elem = max_elem +1;
    min_rate = min(rate_CDK);
    min_elem = find(rate_CDK == min_rate);
    min_elem = min_elem + 1;
    rate_arr(i,:) = [i,min_rate, max_rate];
    figure(1)
    hold on
    %yyaxis right
    %plot(time_vals, cell_length_time);
    %ylabel('Cell Length (um)')
    %yyaxis left
    plot(time_vals, CDK_activity);
     time_idx = [time_vals(min_elem), time_vals(max_elem), max(time_vals)];
%     plot_idx = [CDK_activity(min_elem),CDK_activity(max_elem),CDK_activity(end)];
%     res_show = {strcat('Min:',num2str(min_rate)), strcat('Max:', num2str(max_rate)), num2str(i)};
%     text(time_idx,plot_idx,res_show)
    ylabel('N/C Ratio')
    xlabel('Time(minutes)')
    
    

%hold off
%Protein 2 (Cdc13)


    
    %cell_length_time = (cell2mat(Fluor_cell(i ,2)))*pixel_size;
    
    cell_fluor_time_2 = cell2mat(Fluor_cell(i ,5));
    nuc_fluor_time_2 = cell2mat(Fluor_cell(i,6));
    nuc_fluor_smooth_2 = smoothdata(nuc_fluor_time_2,'movmedian',10);
    time_vals_2 = ((1:length(cell_fluor_time_2))-1)*time_int;
    %CDK_activity = nuc_fluor_time./cell_fluor_time;
    figure(2)
    hold on
    %yyaxis right
    %plot(time_vals, cell_length_time);
    %ylabel('Cell Length (um)')
    %yyaxis left
    plot(time_vals_2, nuc_fluor_smooth_2, time_vals_2, nuc_fluor_time_2);
    text(max(time_vals_2),(nuc_fluor_time_2(end)),num2str(i))
    ylabel(' Mean Cdc 13 Nuclear Intensity (a.u)')
    xlabel ('Time(minutes)')
    
    figure(3)
    hold on
    subplot(1,2,2)
    txt = [num2str(round(div_length(i,3),2))];
    plot(nuc_fluor_smooth_2/max(nuc_fluor_smooth_2), CDK_activity,'DisplayName',txt)
    %plot((nuc_fluor_smooth_2-min(nuc_fluor_smooth_2))/(max(nuc_fluor_smooth_2)-min(nuc_fluor_smooth_2)), CDK_activity,'DisplayName',txt)
    %plot(nuc_fluor_time_2, CDK_activity,'DisplayName',txt)
    xlabel('Mean Nuclear Cdc13 Intensity (a.u.)')
    ylabel('CDK activity')
    title('Smoothed Trace')
    hold on
    subplot(1,2,1)
    %plot(nuc_fluor_time_2/max(nuc_fluor_time_2), CDK_activity,'DisplayName',txt)
    plot(nuc_fluor_time_2, CDK_activity,'DisplayName',txt)
    xlabel('Mean Nuclear Cdc13 Intensity (a.u.)')
    ylabel('CDK activity')
    title('Raw Trace')
    
    %text(max(time_vals_2),(nuc_fluor_time_2(end)),num2str(i))
end
hold off
legend show


%% 
cell_num =2;
wt_max = 2.0;
wt_min = 1.1;
    pixel_size = 0.130;
    time_int = 10; %minutes
    %max_time_frames = 41;
    cell_length_time = (cell2mat(Fluor_cell(cell_num ,2)))*pixel_size;
    fluor_time_1 = cell2mat(Fluor_cell(cell_num ,3));
    fluor_time_2 = cell2mat(Fluor_cell(cell_num ,4));
    
    CDK_activity_time = fluor_time_2./fluor_time_1;
    norm_CDK_Wt = ( CDK_activity_time- wt_min)/(wt_max -wt_min);
    fluor_time_3 = cell2mat(Fluor_cell(cell_num ,6));
    fluor_time_4 = cell2mat(Fluor_cell(cell_num ,5));
    fluor_time_5 = cell2mat(Fluor_cell(cell_num, 7));
    norm_Cdc13_Cyto =  (fluor_time_4-min(fluor_time_4))/(max(fluor_time_4)-min(fluor_time_4));
    fluor_smooth_3 = smoothdata(fluor_time_3,'movmedian',7);
    fluor_smooth_4 = smoothdata(fluor_time_4,'movmedian',7);
    time_vals = ((1:length(fluor_time_3))-1)*time_int;
    %time_vals = ((1:max_time_frames)-1)*time_int;
figure,
    hold on
    %
    plot(time_vals,CDK_activity_time,'Color',[0.6350 0.0780 0.1840],'LineStyle','-','LineWidth',1.0);
   ylabel('Sensor Activity')
   ylim([1.1,2.0])
%     plot(time_vals/max(time_vals),cell_length_time,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0);
%     ylabel('Cell Length')
%     yyaxis left
%     plot(time_vals, fluor_time_3,'Color',[0 0.4470 0.7410],'LineStyle','-','LineWidth',1.0)%
%     plot(time_vals, fluor_time_4,'Color',[0.4660 0.6740 0.1880],'LineStyle','-','LineWidth',1.0)%
%      plot(time_vals, fluor_time_5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0)%
%     ylabel('Cdc13 Mean Intensity')
%     xlabel('Time (minutes)')
     hold off
    
     figure, 
    hold on
    plot(fluor_time_4,  CDK_activity_time)
    %text(fluor_time_3, CDK_activity_time, num2str(cell_length_time));
    xlabel('Mean Cytoplasmic Cdc13 Intensity (a.u.) Normalized')
    ylabel('CDK activity (C/N) Normalized');
    ylim([wt_min, wt_max])
    hold off
    figure, 
    yyaxis left
    plot(time_vals/max(time_vals), fluor_time_4)
     hold on
    yyaxis right
    plot(time_vals/max(time_vals),CDK_activity_time,'Color',[0.6350 0.0780 0.1840],'LineStyle','-','LineWidth',1.0);
   ylabel('Sensor Activity')
   figure, 
      hold on
    yyaxis left
   plot(time_vals,cell_length_time,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0);
    yyaxis right
    plot(time_vals, fluor_time_3,'Color',[0 0.4470 0.7410],'LineStyle','-','LineWidth',1.0)%

    figure, 
      hold on
    yyaxis left
   plot(time_vals,cell_length_time,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0);
    yyaxis right
    plot(time_vals, CDK_activity_time,'Color',[0 0.4470 0.7410],'LineStyle','-','LineWidth',1.0)%
    ylim([1.2,2.8])
    %% 
%      save('Fluorescent_data.mat', 'Fluor_cell')
%  save('Cell_size.mat', 'div_length');
%  save('div_time.mat','div_time');

