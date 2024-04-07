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

max_num_frames = 96;%120 usually
time_int = 5; %minutes
pixel_size = 0.130;
div_time_cut = 1;%usually 1 if you only want cells born in the video
length_ext_cut = 2;
per_pixels = 0.15;
Syn_thresh = 1.10; %  1.14 by default to remove crap traces
time_pts_cut = 10;
time_pts_cut_syn = 11;
time_pts_length = 5;

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
MCM_SYN = zeros(num_cells,17);
figure(1)
figure(2)
figure(3)
for i=1:num_cells
    
    
    cell_length_time = (cell2mat(Fluor_cell(i ,2)))*pixel_size;
    Cyto_Syn = cell2mat(Fluor_cell(i ,3));
    Nuc_Syn = cell2mat(Fluor_cell(i,4));
    
    CDK_activity_Syn_time = Nuc_Syn./ Cyto_Syn ;
    Cdc13_Cyto = cell2mat(Fluor_cell(i,5));
    Cdc13_Nuc = cell2mat(Fluor_cell(i,6));
    Total_cell = cell2mat(Fluor_cell(i, 7));
    CDK_Syn_smooth = smoothdata(CDK_activity_Syn_time,'gaussian',3);
    Cdc13_Cyto_smooth = smoothdata(Cdc13_Cyto,'gaussian',3);
    Cdc13_Nuc_smooth = smoothdata(Cdc13_Nuc, 'gaussian',3);
    Total_cell_smooth = smoothdata(Total_cell, 'gaussian',3);
  
    
    time_vals = ((1:length(cell_length_time ))-1)*time_int;
    [max_Cdc13, id_Cdc13]= max(Cdc13_Nuc_smooth);
    [max_SYN, id_SYN] = max(CDK_Syn_smooth);
    if id_SYN<=time_pts_cut || max_SYN<Syn_thresh || id_SYN + time_pts_length > length(cell_length_time)|| id_SYN<=time_pts_cut_syn
        continue 
    end
    Cdc13_Nuc_Smooth_select =Cdc13_Nuc_smooth(id_SYN - time_pts_cut_syn:id_SYN);
    CDK_SYN_Smooth_select = CDK_Syn_smooth(id_SYN - time_pts_cut_syn:id_SYN);
    cell_length_select = cell_length_time(id_SYN + time_pts_length - time_pts_cut:id_SYN+time_pts_length);
  
    ipt_Cdc13 = findchangepts(Cdc13_Nuc_Smooth_select, 'Statistic', 'linear','MaxNumChanges',1,'MinDistance',5);%2 and 4
    ipt_SYN = findchangepts(CDK_SYN_Smooth_select, 'Statistic', 'linear','MaxNumChanges',1,'MinDistance',5);%2 and 4
    ipt_length = findchangepts(cell_length_select , 'Statistic', 'linear','MaxNumChanges',2,'MinDistance',4);
    loc_Cdc13 = id_SYN - (time_pts_cut_syn - ipt_Cdc13(1) + 1);
    loc_Syn = id_SYN - (time_pts_cut_syn - ipt_SYN(1) + 1);
    loc_length = id_SYN + time_pts_length - (time_pts_cut - ipt_length(1) + 1);
    time_diff = (loc_Syn - loc_Cdc13)*time_int;
    time_diff_growth = (loc_length - loc_Syn)*time_int;
    time_scale = (length(cell_length_time)-1)*time_int;
    MCM_SYN(i,1) = time_diff;
    MCM_SYN (i,2) = ((loc_Cdc13-1)*time_int)/time_scale;% scale so that time zero is beginning of trace (i.e. time since birth)
    MCM_SYN (i,3) = ((loc_Syn-1)*time_int)/time_scale;
    MCM_SYN (i,4) = cell_length_time(loc_Cdc13);
    MCM_SYN (i,5) = cell_length_time(loc_Syn);
    MCM_SYN (i,6) = div_length(i,3) - cell_length_time(loc_Cdc13);
    MCM_SYN (i,7) = div_length(i,3) - cell_length_time(loc_Syn);
    MCM_SYN (i,8) = div_length(i,3);
    MCM_SYN (i,9) = div_length(i,2);
    MCM_SYN (i,10) = MCM_SYN(i,2)* time_scale;
    MCM_SYN (i,11) = MCM_SYN(i,3)*time_scale;
    MCM_SYN (i,12) = time_scale ;
    MCM_SYN (i,13) = Cdc13_Nuc_smooth(loc_Cdc13);
    MCM_SYN (i,14) = Cdc13_Cyto_smooth(loc_Syn);
    MCM_SYN (i,15) = Cdc13_Nuc_smooth (loc_Syn);
    
    
    figure(1)
    hold on

    plot(time_vals, CDK_activity_Syn_time);

    ylabel('N/C Ratio')
    xlabel('Time(minutes)')
    
    
    figure(2)
    hold on
    %yyaxis right
    %plot(time_vals, cell_length_time);
    %ylabel('Cell Length (um)')
    %yyaxis left
    plot(time_vals, Cdc13_Nuc);
    text(max(time_vals),(Cdc13_Nuc(end)),num2str(i))
    ylabel(' Mean Cdc 13 Nuclear Intensity (a.u)')
    xlabel ('Time(minutes)')
    
    figure(3)
    hold on
    subplot(1,2,2)
    txt = [num2str(round(div_length(i,3),2))];
    plot((Cdc13_Cyto-min(Cdc13_Cyto))/(max(Cdc13_Cyto)-min(Cdc13_Cyto)), (CDK_activity_Syn_time - min(CDK_activity_Syn_time))/(max(CDK_activity_Syn_time)-min(CDK_activity_Syn_time)),'DisplayName',txt)
    
    xlabel('Normalized mean Cdc13 Cytoplasmic Intensity (a.u.)')
    ylabel('N/C Ratio Normalized')
    title(' Normalized')
    hold on
    subplot(1,2,1)
    %plot(nuc_fluor_time_2/max(nuc_fluor_time_2), CDK_activity,'DisplayName',txt)
    plot(Cdc13_Cyto, CDK_activity_Syn_time,'DisplayName',txt)
    xlabel('Mean Cytoplasmic Cdc13 Intensity (a.u.)')
    ylabel('N/C Ratio')
    title('Raw Trace')
    
    %text(max(time_vals_2),(nuc_fluor_time_2(end)),num2str(i))
end
hold off
legend show
elem_nonzeros = find(MCM_SYN(:,2));
MCM_SYN = MCM_SYN(elem_nonzeros,:);
Fluor_cell = Fluor_cell(elem_nonzeros,:);
div_length = div_length(elem_nonzeros,:);
div_time = div_time(elem_nonzeros,:);
cov_Cdc13 = std(MCM_SYN(:,4))/mean(MCM_SYN(:,4));
cov_SYN = std(MCM_SYN(:,5))/mean(MCM_SYN(:,5));
MCM_SYN(1,16) = cov_Cdc13;
MCM_SYN(1,17) = cov_SYN;

%% 
cell_num =4;
time_pts_cut = 10;
time_pts_cut_syn = 11;
time_pts_length = 5;

    cell_length_time = (cell2mat(Fluor_cell(cell_num ,2)))*pixel_size;
    Cyto_Syn = cell2mat(Fluor_cell(cell_num ,3));
    Nuc_Syn = cell2mat(Fluor_cell(cell_num ,4));
    
    CDK_activity_Syn_time = Nuc_Syn./ Cyto_Syn ;
    Cdc13_Cyto = cell2mat(Fluor_cell(cell_num ,5));
    Cdc13_Nuc = cell2mat(Fluor_cell(cell_num ,6));
    Total_cell = cell2mat(Fluor_cell(cell_num, 7));
    CDK_Syn_smooth = smoothdata(CDK_activity_Syn_time,'gaussian',3);
    Cdc13_Cyto_smooth = smoothdata(Cdc13_Cyto,'gaussian',3);
    Cdc13_Nuc_smooth = smoothdata(Cdc13_Nuc, 'gaussian',3);
    Total_cell_smooth = smoothdata(Total_cell, 'gaussian',3);
    
    time_vals = ((1:length(cell_length_time ))-1)*time_int;
    [max_Cdc13, id_Cdc13]= max(Cdc13_Nuc_smooth);
    [max_SYN, id_SYN] = max(CDK_Syn_smooth);
    Cdc13_Nuc_Smooth_select =Cdc13_Nuc_smooth(id_SYN - time_pts_cut_syn:id_SYN);
    CDK_SYN_Smooth_select = CDK_Syn_smooth(id_SYN - time_pts_cut_syn:id_SYN);
    cell_length_select = cell_length_time(id_SYN + time_pts_length - time_pts_cut:id_SYN+time_pts_length);

    ipt_Cdc13 = findchangepts(Cdc13_Nuc_Smooth_select, 'Statistic', 'linear','MaxNumChanges',1,'MinDistance',5);%2 and 4
    ipt_SYN = findchangepts(CDK_SYN_Smooth_select, 'Statistic', 'linear','MaxNumChanges',1,'MinDistance',5);%2 and 4
    ipt_length = findchangepts(cell_length_select , 'Statistic', 'linear','MaxNumChanges',2,'MinDistance',4);
    loc_Cdc13 = id_SYN - (time_pts_cut_syn - ipt_Cdc13(1) + 1);
    loc_Syn = id_SYN - (time_pts_cut_syn - ipt_SYN(1) + 1);
    loc_length = id_SYN + time_pts_length - (time_pts_cut - ipt_length(1) + 1);
    figure,
    hold on
    yyaxis left
%     plot(time_vals/max(time_vals),CDK_activity_time_syn/max(CDK_activity_time_syn),'b-', time_vals/max(time_vals), CDK_activity_time_mcm/max(CDK_activity_time_mcm),'r-');
    plot(time_vals/max(time_vals),Cdc13_Nuc_smooth,'b-',time_vals/max(time_vals),Cdc13_Nuc,'b--','LineWidth',1)
    text(time_vals(loc_Cdc13)/max(time_vals), Cdc13_Nuc_smooth(loc_Cdc13),'Nuclear Export')
    hold on
    ylabel('Nuclear Mean Intensity')
     yyaxis right
    plot(time_vals/max(time_vals),CDK_Syn_smooth,'r-',time_vals/max(time_vals),CDK_activity_Syn_time,'r--','LineWidth',1)
    text(time_vals(loc_Syn)/max(time_vals), CDK_Syn_smooth(loc_Syn),'Cytoplasmic Entry')
   % plot(time_vals/max(time_vals), (CDK_activity_time_mcm - min(CDK_activity_time_mcm))/(max(CDK_activity_time_mcm)- min(CDK_activity_time_mcm)))
    ylabel('Sensor Activity')
    hold off
    time_vals_Cdc13 = time_vals(id_SYN - time_pts_cut_syn:id_SYN);
    time_vals_syn = time_vals(id_SYN - time_pts_cut_syn:id_SYN);
    time_vals_length = time_vals(id_SYN+time_pts_length - time_pts_cut:id_SYN+time_pts_length);
    figure, 
    plot (time_vals_syn,CDK_SYN_Smooth_select);%,'b-', time_vals, CDK_activity_mcm_smooth,'g-')
    text(time_vals_syn(ipt_SYN), CDK_SYN_Smooth_select(ipt_SYN),'Rate Change SYN')
    figure
    plot (time_vals_Cdc13,Cdc13_Nuc_Smooth_select);
    text(time_vals_Cdc13(ipt_Cdc13),Cdc13_Nuc_Smooth_select(ipt_Cdc13),'Rate Change Cdc13')
    figure, 
    plot (time_vals_length,cell_length_select);
    text(time_vals_length(ipt_length),cell_length_select(ipt_length),'Rate Change Cell Length')
%% 
figure,
    hold on
    yyaxis right
    plot(time_vals/max(time_vals),CDK_Syn_smooth,'r-',time_vals/max(time_vals), CDK_activity_Syn_time,'r--','LineWidth',1.0);
   ylabel('Sensor Activity')
   
%     plot(time_vals/max(time_vals),cell_length_time,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0);
%     ylabel('Cell Length')
    yyaxis left
    plot(time_vals/max(time_vals), Cdc13_Nuc_smooth,'b-',time_vals/max(time_vals), Cdc13_Nuc_smooth,'b--','LineWidth',1.0)%
    plot(time_vals/max(time_vals), Cdc13_Cyto_smooth,'m-',time_vals/max(time_vals), Cdc13_Cyto_smooth,'m--','LineWidth',1.0)%
    plot(time_vals/max(time_vals), Total_cell_smooth,'g-',time_vals/max(time_vals), Total_cell_smooth,'g--','LineWidth',1.0)%
    ylabel('Cdc13 Mean Intensity')
    xlabel('Time (minutes)')
    hold off
    
     figure, 
    hold on
    norm_Cdc13_cyto  = (Cdc13_Cyto -min(Cdc13_Cyto))/(max(Cdc13_Cyto)-min(Cdc13_Cyto));
    norm_Cdc13_Cyto_Smooth = (Cdc13_Cyto_smooth -min(Cdc13_Cyto_smooth))/(max(Cdc13_Cyto_smooth)-min(Cdc13_Cyto_smooth));
    norm_Syn_smooth = (CDK_Syn_smooth - min(CDK_Syn_smooth))/(max(CDK_Syn_smooth)-min(CDK_Syn_smooth));
    norm_Syn = (CDK_activity_Syn_time - min( CDK_activity_Syn_time))/(max(CDK_activity_Syn_time)-min(CDK_activity_Syn_time));
    plot(norm_Cdc13_cyto, norm_Syn)
    %(CDK_activity_syn_smooth - min(CDK_activity_syn_smooth))/(max(CDK_activity_syn_smooth)
    %text(fluor_time_3, CDK_activity_time, num2str(cell_length_time));
    xlabel('Mean Cytoplasmic Cdc13 Intensity (a.u.)')
    ylabel('N/C Ratio');
     title('Phase Plot')
    hold off
    figure, 
    yyaxis left
    plot(time_vals/max(time_vals), Cdc13_Cyto_smooth,'m-',time_vals/max(time_vals), Cdc13_Cyto,'m--','LineWidth',1.0)%
     hold on
    yyaxis right
    plot(time_vals/max(time_vals),CDK_Syn_smooth,'r-',time_vals/max(time_vals), CDK_activity_Syn_time,'r--','LineWidth',1.0);
   ylabel('Sensor Activity')
   figure, 
   yyaxis left
   plot (time_vals, cell_length_time);
   hold on
   yyaxis right
   plot(time_vals, CDK_Syn_smooth)
   figure (10)
   plot(Cdc13_Cyto,norm_Syn);
    hold on
    %% 
%      save('Fluorescent_data.mat', 'Fluor_cell')
%  save('Cell_size.mat', 'div_length');
%  save('div_time.mat','div_time');

