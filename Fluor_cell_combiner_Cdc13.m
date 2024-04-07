% Fluor_cell_combined = [Fluor_cell_01;Fluor_cell_02];
% uniq_cell_arr_combined = [uniq_cell_arr_01;uniq_cell_arr_02];
Fluor_cell = Fluor_cell_combined;
uniq_cell_arr = uniq_cell_arr_combined;
pixel_size = 0.130;
per_pixels = 0.15;
time_pts_cut = 10;
time_pts_length = 5;
time_pts_cut_MCM = 15;
time_int = 5; %minutes
num_cells = length(Fluor_cell(:,1));
MCM_SYN = zeros(num_cells,17);
figure(1)
figure(2)
figure(3)
figure (4)
figure (5)
figure (6)
figure (7)
figure (8)
figure (9)
for i=1:num_cells
   cell_dat = uniq_cell_arr{i};
   cell_length_time = (cell2mat(Fluor_cell(i ,2)))*pixel_size;
   Cdc13_Cyto = cell2mat(cell_dat(:,4));
    Cdc13_Nuc = cell2mat(cell_dat(:,5));
    Cdc13_total_cell = cell2mat(Fluor_cell(i ,3));
    %CDK_activity_MCM_time = Cyto_MCM./ Nuc_MCM ;
%     Cdc13_Cyto = cell2mat(Fluor_cell(i,5));
%     Cdc13_Nuc = cell2mat(Fluor_cell(i,6));
%     Total_cell = cell2mat(Fluor_cell(i, 7));
     Cdc13_Cyto_Smooth = smoothdata(Cdc13_Cyto,'sgolay',3);
    Cdc13_Nuc_Smooth = smoothdata(Cdc13_Nuc ,'sgolay',3);
   Cdc13_total_cell_smooth = smoothdata(Cdc13_total_cell,'sgolay',3);
%     Cdc13_Cyto_smooth = smoothdata(Cdc13_Cyto,'sgolay',4);
%     Cdc13_Nuc_smooth = smoothdata(Cdc13_Nuc,'sgolay',4);
%     Total_cell_smooth = smoothdata(Total_cell, 'sgolay',4);
%     norm_Cdc13_Nuc_raw = (Cdc13_Nuc-min(Cdc13_Nuc))/(max(Cdc13_Nuc)-min(Cdc13_Nuc));
%     norm_MCM_raw = (CDK_activity_MCM_time-min(CDK_activity_MCM_time))/(max(CDK_activity_MCM_time)-min(CDK_activity_MCM_time));
    
    time_vals = ((1:length(Cdc13_Nuc ))-1)*time_int;
    norm_time = time_vals/max(time_vals);
    
    figure(1)
    hold on
    plot(norm_time, Cdc13_Nuc);
    ylabel('Nuclear Intensity')
    xlabel('Cell Cycle Time Normalized')
    figure(2)
    hold on
    plot(norm_time, Cdc13_Cyto);
    ylabel('Cytoplasmic Intensity')
    xlabel('Cell Cycle Time Normalized')
    
    figure (3)
    
    hold on
    plot(norm_time, Cdc13_total_cell);
    ylabel('Total Mean Intensity')
    xlabel('Cell Cycle Time Normalized')

    figure (4)
    hold on
    plot(norm_time, Cdc13_Nuc_Smooth);
    ylabel('Smoothed Nuclear Mean Intensity')
    xlabel('Cell Cycle Time Normalized')
    figure (5)
    hold on
    plot(norm_time, Cdc13_Cyto_Smooth);
    ylabel('Smoothed Cytoplasmic Mean')
    xlabel('Cell Cycle Time Normalized')
    figure (6)
    hold on
    plot(norm_time, Cdc13_total_cell_smooth);
    ylabel('Smoothed Nuclear Intensity')
    xlabel('Cell Cycle Time Normalized')
    figure (7)
     subplot(3,1,2)
    plot(time_vals/max(time_vals), Cdc13_Nuc_Smooth,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0);
    title('Nucleus')
    hold on
    subplot(3,1,1)
    plot(time_vals/max(time_vals),  Cdc13_total_cell_smooth,'Color',[0.4660 0.6740 0.1880],'LineStyle','-','LineWidth',1.0);
    title('Total Mean Cell Intensity')
    hold on
%     subplot(4,1,2)
%     plot(time_vals/max(time_vals), cell_tot_abs_time,'Color',[0.6350 0.0780 0.1840],'LineStyle','-','LineWidth',1.0);
%     title('Total Absolute Intensity')
%     hold on
    subplot(3,1,3)
    plot(time_vals/max(time_vals), Cdc13_Cyto_Smooth , 'Color',[0 0.4470 0.7410],'LineStyle','-','LineWidth',1.0);
    title('Cytoplasm')
    hold on
   
%     figure(2)
%     hold on
% 
%     plot(time_vals, Cdc13_Nuc);
%     text(max(time_vals),(Cdc13_Nuc(end)),num2str(i))
%     ylabel(' Mean Cdc 13 Nuclear Intensity (a.u)')
%     xlabel ('Time(minutes)')
%     
%     figure(3)
%     hold on
%     subplot(1,2,2)
%     norm_Cdc13_nuc = (Cdc13_Nuc_smooth-min(Cdc13_Nuc_smooth))/(max(Cdc13_Nuc_smooth)-min(Cdc13_Nuc_smooth));
%     norm_MCM = (CDK_MCM_smooth-min(CDK_MCM_smooth))/(max(CDK_MCM_smooth)-min(CDK_MCM_smooth));
%     plot(norm_Cdc13_nuc,norm_MCM)
%     xlabel('Mean Normalized Nuclear Cdc13 Intensity (a.u.)')
%     ylabel('C/N Ratio Normalized')
%     title('Normalized Smoothed Trace')
%     hold on
%     subplot(1,2,1)
%     plot(Cdc13_Nuc, CDK_activity_MCM_time)
%     xlabel('Mean Nuclear Cdc13 Intensity (a.u.)')
%     ylabel('C/N Ratio')
%     title('Raw Trace')
%     figure (4)
%     hold on
%     plot(norm_Cdc13_nuc,norm_MCM)
     
%     figure (6)
%     hold on
%     plot(norm_time, Cdc13_Nuc)
%     figure (7)
%     hold on
%     plot(norm_time, Cdc13_Cyto);
%     figure (8)
%     hold on
%     plot(norm_time, Total_cell)
%     figure (9)
%     hold on
% 
%     plot(Cdc13_Nuc_smooth, norm_MCM)
%     %text(max(time_vals_2),(nuc_fluor_time_2(end)),num2str(i))
end
hold off
%% 


time_int = 5;
cell_num =  12;
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
    Cdc13_Cyto_smooth = smoothdata(cell_fluor_time,'sgolay',3);
    Cdc13_Nuc_smooth = smoothdata(nuc_fluor_time ,'sgolay',3);
   Cdc13_total_cell_smooth = smoothdata(cell_mean,'sgolay',3);
%     cell_fluor_smooth = smoothdata(cell_fluor_time,'movmedian',med_filter_sz);
%     nuc_fluor_smooth = smoothdata(nuc_fluor_time,'movmedian',med_filter_sz);
    time_vals = time_points*time_int;
%     [max_elem,Imax] = max(cell_fluor_time(15:end));
%     Imax_min = Imax - 1;
%     nuc_diam_max = nuc_diam_time(Imax);
%     nuc_diam_min = nuc_diam_time (Imax_min);
%     nuc_fluor_max = nuc_fluor_time(Imax);
%     nuc_fluor_min = nuc_fluor_time(Imax_min);
%     cell_fluor_max = cell_fluor_time(Imax);
%     cell_fluor_min = cell_fluor_time(Imax_min);
%     cell_length_max = cell_length_time(Imax);
%     cell_length_min = cell_length_time(Imax_min);
%     nuc_bg_intensities_mean = min(nuc_fluor_time);
%      cell_bg_intensitie_mean = min(cell_fluor_time);
%       nuc_abs_min = (4/3)*pi*((nuc_diam_min/2)^3) * (nuc_fluor_min - nuc_bg_intensities_mean);
%     nuc_abs_max = (4/3)*pi*((nuc_diam_max/2)^3) * (nuc_fluor_max - nuc_bg_intensities_mean);
%     cell_abs_min = ((3.5^2)* cell_length_min*pi - (4/3)*pi*((nuc_diam_min/2)^3))*(cell_fluor_min-cell_bg_intensitie_mean);
%     cell_abs_max = ((3.5^2)* cell_length_max*pi - (4/3)*pi*((nuc_diam_max/2)^3))*(cell_fluor_max- cell_bg_intensitie_mean);
%     nuc_abs_min_alt = ((3.5/2)^2)* cell_length_min*pi * 0.08 *nuc_fluor_min;
%     nuc_abs_max_alt = ((3.5/2)^2)* cell_length_max*pi * 0.08 *nuc_fluor_max;
%     diff_nuc  = nuc_abs_max -nuc_abs_min;
%     diff_nuc_alt = nuc_abs_max_alt - nuc_abs_min_alt;
%     diff_cyto = cell_abs_max - cell_abs_min;
%     conc_nuc_diff = nuc_fluor_max - nuc_fluor_min;
%     conc_cyto_diff = cell_fluor_max - cell_fluor_min;
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
%    figure, 
%    scatter(diff_array_conc(:,1), diff_array_conc(:,2));
   figure, 
   plot(time_vals/max(time_vals), cell_fluor_time./nuc_fluor_time,'b-');
   figure, 
   subplot(3,1,3)
   plot(time_vals/max(time_vals), cell_fluor_time, 'Color',[0 0.4470 0.7410],'LineStyle','-','LineWidth',1.0);
    title('Cytoplasmic')
    ylabel('Mean Intensity (A.U)')
    xlabel('Time(normalized)')
%     xlim([0.6,0.80])
%      ylim([230,270])
   hold on 
   subplot(3,1,2)
   plot(time_vals/max(time_vals), nuc_fluor_time, 'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0);
    title('Nuclear')
    ylabel('Mean Intensity (A.U)')
     xlabel('Time(normalized)')
%      xlim([0.6,0.80])
%      ylim([300,700])
    hold on
    subplot(3,1,1);
    plot(time_vals/max(time_vals), cell_mean, 'Color',[0.4660 0.6740 0.1880],'LineStyle','-','LineWidth',1.0);
    title('Total Cell')
    ylabel('Mean Intensity (A.U)')
     xlabel('Time(normalized)')
%      xlim([0.6,0.80])
%      ylim([225,300])
     hold off
%    figure, 
%    subplot(3,1,3)
%    plot(time_vals/max(time_vals),  Cdc13_Cyto_smooth, 'Color',[0 0.4470 0.7410],'LineStyle','-','LineWidth',1.0);
%     title('Cytoplasmic')
%     ylabel('Mean Intensity (A.U)')
%     xlabel('Time(normalized)')
%    hold on 
%    subplot(3,1,2)
%    plot(time_vals/max(time_vals),Cdc13_Nuc_smooth, 'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0);
%     title('Nuclear Smooth')
%     ylabel('Mean Intensity (A.U)')
%      xlabel('Time(normalized)')
%     hold on
%     subplot(3,1,1);
%     plot(time_vals/max(time_vals), Cdc13_total_cell_smooth, 'Color',[0.4660 0.6740 0.1880],'LineStyle','-','LineWidth',1.0);
%     title('Total Cell Smoothed')
%     ylabel('Mean Intensity (A.U)')
%      xlabel('Time(normalized)')