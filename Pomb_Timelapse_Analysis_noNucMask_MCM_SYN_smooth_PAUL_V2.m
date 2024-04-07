%% 
[filenames_fluor,path_fluorescence] = uigetfile ('.tif','Pick SynCut Fluorescent Files','Multiselect','on');
[filenames_Prot2_fluor,path_Prot2_fluor] = uigetfile ('.tif','Pick MCM Fluorescent Files','Multiselect','on');
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

max_num_frames =50;%120 usually
time_int =5; %minutes
pixel_size = 0.130;
div_time_cut = 1;%usually 1 if you only want cells born in the video
length_ext_cut = 3;
per_pixels = 0.15;
time_pts_cut_mcm = 15;%default 13
time_pts_cut_syn = 10;%default 10
time_pts_length = 5;%default 7
time_pts_length_cut = 10;%default 10
Syn_thresh = 1.14; % default 1.14 in case you have a mixed population where some cells don't have syncut3.

BD_data = BD_file.data;
elem_born = find(BD_data(:,2)>div_time_cut);
BD_data_revised = BD_data(elem_born,:);
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
% find_div_errs = find(div_time(:,2)>1);
% div_time = div_time(find_div_errs,:);
% div_length = div_length(find_div_errs,:);
length_extension = div_length(:,3)-div_length(:,2);
length_err = find(length_extension > length_ext_cut);
 div_time = div_time(length_err,:);
 div_length = div_length(length_err,:);
 length_extension = length_extension(length_err,:);

figure (4),

histogram(div_time(:,2));
mean(div_time(:,2));
figure(6)
histogram(div_length(:,3));
mean_div = mean(div_length(:,3));
std_div = std(div_length(:,3));
cov_value = std_div/mean_div
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
Fluor_cell = cell(num_cells, 6); %cell_ID, cell_length, cell intensity, nuclear intensity
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
    for j = cell_born_frame:cell_death_frame
      frame_grab_Prot2 = imread(strcat(path_Prot2_fluor,filenames_Prot2_fluor{j}));
      frame_grab_fluor = imread(strcat(path_fluorescence,filenames_fluor{j}));
      %frame_grab_bin = imread(strcat(path_bin,filenames_bin{j}));
      frame_grab_trk_bin = imread(trk_image_files{j});
%       nuc_stats = regionprops(logical(frame_grab_nuc),'Centroid');
%       nuc_centroids = cell2mat({nuc_stats.Centroid}');
      cell_stats = regionprops(logical(frame_grab_trk_bin),frame_grab_trk_bin, 'PixelList','MeanIntensity','MaxFeretProperties');
      cell_stats_fluor = regionprops(logical(frame_grab_trk_bin),frame_grab_fluor, 'PixelValues','PixelList','MeanIntensity','MaxFeretProperties');
      cell_stats_fluor_2 = regionprops(logical(frame_grab_trk_bin),frame_grab_Prot2, 'PixelValues','PixelList','MeanIntensity','MaxFeretProperties');
    
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
    
end
Fluor_cell = Fluor_cell(~cellfun(@isempty, Fluor_cell(:,1)), :);
num_cells = length(Fluor_cell(:,1));
%% Protein 1(Syncut)
rate_arr = zeros(num_cells,3);
% figure(1)
% figure(2)
% figure(3)
MCM_SYN = zeros(num_cells,17);
figure(8)
for i=1:num_cells
    
    cell_length_time = (cell2mat(Fluor_cell(i ,2)))*pixel_size;
    
    cell_fluor_time = cell2mat(Fluor_cell(i ,3));
    nuc_fluor_time = cell2mat(Fluor_cell(i,4));
    cell_fluor_time_2 = cell2mat(Fluor_cell(i ,5));
    nuc_fluor_time_2 = cell2mat(Fluor_cell(i,6));
    %nuc_fluor_smooth = smoothdata(nuc_fluor_time,'movmedian',3);
    time_vals = ((1:length(cell_fluor_time))-1)*time_int;
    CDK_activity_syn = nuc_fluor_time./cell_fluor_time;
    CDK_activity_mcm = cell_fluor_time_2./nuc_fluor_time_2;
    CDK_activity_syn_smooth = smoothdata(CDK_activity_syn, 'sgolay',4);
  %CDK_activity_syn_smooth = smoothdata(CDK_activity_time_syn, 'gaussian',3);
    CDK_activity_mcm_smooth = smoothdata(CDK_activity_mcm, 'sgolay',4);
    %fluor_time_3 = cell2mat(Fluor_cell(cell_num ,6));
    %fluor_smooth_3 = smoothdata(fluor_time_3,'movmedian',5);
    %time_vals = ((1:length(fluor_time_2))-1)*time_int;
    [max_MCM, id_MCM]= max(CDK_activity_mcm_smooth);
    [max_SYN, id_SYN] = max(CDK_activity_syn_smooth);
    if id_SYN<=time_pts_cut_syn || max_SYN<Syn_thresh|| id_MCM + time_pts_length > length(cell_length_time)||id_MCM<=time_pts_cut_mcm
        continue 
    end
    CDK_MCM_Smooth_select = CDK_activity_mcm_smooth(id_MCM - time_pts_cut_mcm:id_MCM);
    CDK_SYN_Smooth_select = CDK_activity_syn_smooth(id_SYN - time_pts_cut_syn:id_SYN);
    
   
     cell_length_select = cell_length_time(id_MCM+time_pts_length - time_pts_length_cut:id_MCM+time_pts_length);
    %ipt_MCM = findchangepts(CDK_MCM_Smooth_select , 'Statistic', 'linear','MinThreshold',0.00065);%1 and 4
    ipt_MCM = findchangepts(CDK_MCM_Smooth_select , 'Statistic', 'linear','MinThreshold',0.002,'MinDistance',4);%1 and 4
    ipt_SYN = findchangepts(CDK_SYN_Smooth_select, 'Statistic', 'linear','MinThreshold',0.001,'MinDistance',4);%1 and 3
    ipt_length = findchangepts(cell_length_select , 'Statistic', 'linear','MaxNumChanges',2,'MinDistance',4);% 1 and 4
    loc_MCM = id_MCM - (time_pts_cut_mcm - ipt_MCM(1) + 1);
    loc_Syn = id_SYN - (time_pts_cut_syn - ipt_SYN(1) + 1);
    loc_length = id_MCM+time_pts_length - (time_pts_length_cut  - ipt_length(1) + 1);
    time_diff = (loc_Syn - loc_MCM)*time_int;
    time_diff_growth = (loc_length - loc_MCM)*time_int;
    time_scale = (length(cell_length_time)-1)*time_int;
    MCM_SYN(i,1) = time_diff;
    MCM_SYN (i,2) = ((loc_MCM-1)*time_int)/time_scale;% scale so that time zero is beginning of trace (i.e. time since birth)
    MCM_SYN (i,3) = ((loc_Syn-1)*time_int)/time_scale;
    MCM_SYN (i,4) = cell_length_time(loc_MCM);
    MCM_SYN (i,5) = cell_length_time(loc_Syn);
    MCM_SYN (i,6) = div_length(i,3) - cell_length_time(loc_MCM);
    MCM_SYN (i,7) = div_length(i,3) - cell_length_time(loc_Syn);
    MCM_SYN (i,8) = div_length(i,3);
    MCM_SYN (i,9) = div_length(i,2);
    MCM_SYN (i,10) = MCM_SYN(i,2)* time_scale;
    MCM_SYN (i,11) = MCM_SYN(i,3)*time_scale;
    MCM_SYN (i,12) = time_scale ;
    MCM_SYN (i,13) = CDK_activity_mcm_smooth(loc_MCM);
    MCM_SYN (i,14) = time_diff_growth;
    MCM_SYN (i,15) = CDK_activity_mcm_smooth (loc_Syn);
    %MCM_SYN = MCM_SYN(MCM_SYN>0);
    figure (8)
    plot(time_vals/max(time_vals),CDK_activity_syn);
    hold on

end
hold off
% legend show
elem_nonzeros = find(MCM_SYN(:,2));
MCM_SYN = MCM_SYN(elem_nonzeros,:);
Fluor_cell = Fluor_cell(elem_nonzeros,:);
cov_MCM = std(MCM_SYN(:,4))/mean(MCM_SYN(:,4));
cov_SYN = std(MCM_SYN(:,5))/mean(MCM_SYN(:,5));
MCM_SYN(1,16) = cov_MCM;
MCM_SYN(1,17) = cov_SYN;

MCM_change = mean(MCM_SYN(:,2));
SYN_change = mean(MCM_SYN(:,3));
[b2,Sfit_MCM] = polyfit(MCM_SYN(:,4), MCM_SYN(:,8),1);
[Yfit_MCM, delta_fit_MCM] = polyconf(b2, MCM_SYN(:,4), Sfit_MCM);
mdl_MCM =  fitlm(MCM_SYN(:,4),MCM_SYN(:,8));
yCalc1_MCM = polyval(b2, MCM_SYN(:,4));
figure(5), 
scatter(MCM_SYN(:,4), MCM_SYN(:,8))
hold on 
plot(MCM_SYN(:,4), yCalc1_MCM );
hold off
%% 

cell_num=37;



time_pts_cut_mcm = 15;%default 13
time_pts_cut_syn = 10;%default 10
time_pts_length = 5;%default 7
time_pts_length_cut = 10;%default 10
    cell_length_time = (cell2mat(Fluor_cell(cell_num ,2)))*pixel_size;
    fluor_time_1 = cell2mat(Fluor_cell(cell_num ,3));
    fluor_time_2 = cell2mat(Fluor_cell(cell_num ,4));
    fluor_time_3 = cell2mat(Fluor_cell(cell_num ,5));
    fluor_time_4 = cell2mat(Fluor_cell(cell_num ,6));
    CDK_activity_time_syn = fluor_time_2./fluor_time_1;
    CDK_activity_time_mcm = fluor_time_3./fluor_time_4;
    CDK_activity_syn_smooth = smoothdata(CDK_activity_time_syn, 'sgolay',4);
  %CDK_activity_syn_smooth = smoothdata(CDK_activity_time_syn, 'gaussian',3);
    CDK_activity_mcm_smooth = smoothdata(CDK_activity_time_mcm, 'sgolay',4);
    %fluor_time_3 = cell2mat(Fluor_cell(cell_num ,6));
    %fluor_smooth_3 = smoothdata(fluor_time_3,'movmedian',5);
    time_vals = ((1:length(fluor_time_3))-1)*time_int;
    [max_MCM, id_MCM]= max(CDK_activity_mcm_smooth);
    [max_SYN, id_SYN] = max(CDK_activity_syn_smooth);
    CDK_MCM_Smooth_select = CDK_activity_mcm_smooth(id_MCM - time_pts_cut_mcm :id_MCM);
    CDK_SYN_Smooth_select = CDK_activity_syn_smooth(id_SYN - time_pts_cut_syn:id_SYN);
    cell_length_select = cell_length_time(id_MCM+time_pts_length - time_pts_length_cut:id_MCM+time_pts_length);
%     ipt_MCM = findchangepts(CDK_MCM_Smooth_select , 'Statistic', 'linear','MaxNumChanges',2,'MinDistance',5);%1 and 4
%     ipt_MCM = findchangepts(CDK_MCM_Smooth_select , 'Statistic', 'linear','MinThreshold',0.00065);%1 and 4
     ipt_MCM = findchangepts(CDK_MCM_Smooth_select , 'Statistic', 'linear','MinThreshold',0.002,'MinDistance',4);%1 and 4
    ipt_SYN = findchangepts(CDK_SYN_Smooth_select, 'Statistic', 'linear','MinThreshold',0.001,'MinDistance',4);%1 and 3
    ipt_length = findchangepts(cell_length_select , 'Statistic', 'linear','MaxNumChanges',2,'MinDistance',4);% 1 and 4
    loc_MCM = id_MCM - (time_pts_cut_mcm - ipt_MCM(1) + 1);
    loc_Syn = id_SYN - (time_pts_cut_syn - ipt_SYN(1) + 1);
    loc_length = id_MCM+time_pts_length - (time_pts_length_cut  - ipt_length(1) + 1);
    figure,
    hold on
    yyaxis right
     %plot(time_vals/max(time_vals),CDK_activity_time_syn/max(CDK_activity_time_syn),'b-', time_vals/max(time_vals), CDK_activity_time_mcm/max(CDK_activity_time_mcm),'r-');
    plot(time_vals/max(time_vals),(CDK_activity_syn_smooth - min(CDK_activity_syn_smooth))/(max(CDK_activity_syn_smooth) - min(CDK_activity_syn_smooth)),'r-',time_vals/max(time_vals),(CDK_activity_time_syn - min(CDK_activity_time_syn))/(max(CDK_activity_time_syn) - min(CDK_activity_time_syn)),'r--','LineWidth',1)
    text(time_vals(loc_Syn)/max(time_vals), (CDK_activity_time_syn(loc_Syn) - min(CDK_activity_time_syn))/(max(CDK_activity_time_syn) - min(CDK_activity_time_syn)),'Cytoplasmic Entry')
%   
    hold on
     plot( time_vals/max(time_vals), (CDK_activity_mcm_smooth - min(CDK_activity_mcm_smooth))/(max(CDK_activity_mcm_smooth)- min(CDK_activity_mcm_smooth)),'b-',time_vals/max(time_vals), (CDK_activity_time_mcm- min(CDK_activity_time_mcm))/(max(CDK_activity_time_mcm)- min(CDK_activity_time_mcm)),'b--','LineWidth',1);
     text(time_vals(loc_MCM)/max(time_vals), (CDK_activity_time_mcm(loc_MCM) - min(CDK_activity_time_mcm))/(max(CDK_activity_time_mcm)- min(CDK_activity_time_mcm)),'Nuclear Entry')
    ylabel('Sensor Activity')
     yyaxis left
    plot(time_vals/max(time_vals), cell_length_time,'b-','LineWidth',1);
   % plot(time_vals/max(time_vals), (CDK_activity_time_mcm - min(CDK_activity_time_mcm))/(max(CDK_activity_time_mcm)- min(CDK_activity_time_mcm)))
    ylabel('Cell Length')
    hold off
    time_vals_mcm = time_vals(id_MCM - time_pts_cut_mcm:id_MCM);
    time_vals_syn = time_vals(id_SYN - time_pts_cut_syn:id_SYN);
    time_vals_length = time_vals(id_MCM+time_pts_length - time_pts_length_cut :id_MCM+time_pts_length);
    figure, 
    plot (time_vals_syn,CDK_SYN_Smooth_select);%,'b-', time_vals, CDK_activity_mcm_smooth,'g-')
    text(time_vals_syn(ipt_SYN), CDK_SYN_Smooth_select(ipt_SYN),'Max Rate Change SYN')
    figure
    plot (time_vals_mcm,CDK_MCM_Smooth_select);
    text(time_vals_mcm(ipt_MCM),CDK_MCM_Smooth_select(ipt_MCM),'Max Rate Change MCM')
    figure, 
    plot (time_vals_length,cell_length_select);
    text(time_vals_length(ipt_length),cell_length_select(ipt_length),'Max Rate Change Cell Length')
      figure, 
      plot(time_vals/max(time_vals),CDK_activity_time_syn )
%     hold on
%     yyaxis right
%     plot(time_vals/max(time_vals),CDK_activity_syn_smooth/max(CDK_activity_syn_smooth),'b-', time_vals/max(time_vals), CDK_activity_mcm_smooth/max(CDK_activity_mcm_smooth),'r-');
%     ylabel('Sensor Activity')
%     yyaxis left
%     plot(time_vals/max(time_vals), cell_length_time);
%     ylabel('Cell Length')
%     hold off  
    
    
%     figure, 
%     hold on
%     plot(fluor_time_3, CDK_activity_time)
%     %text(fluor_time_3, CDK_activity_time, num2str(cell_length_time));
%     xlabel('Mean Cdc13 Intensity (a.u.)')
%     ylabel('CDK activity (N/C)');
%     
    hold off

    %% 
%      save('Fluorescent_data.mat', 'Fluor_cell')
%  save('Cell_size.mat', 'div_length');
%  save('div_time.mat','div_time');
%mean(MCM_SYN)
extension = MCM_Combined(:,8)- MCM_Combined(:,9);
[b2,Sfit_MCM] = polyfit(MCM_Combined(:,9), extension,1);
[Yfit_MCM, delta_fit_MCM] = polyconf(b2, extension, Sfit_MCM);
mdl_MCM =  fitlm(MCM_Combined(:,9), extension,'linear');
yCalc1_MCM = polyval(b2, MCM_Combined(:,9));
figure(5), 
scatter(MCM_Combined(:,9),extension)
hold on 
plot(MCM_Combined(:,9), yCalc1_MCM );
hold off
