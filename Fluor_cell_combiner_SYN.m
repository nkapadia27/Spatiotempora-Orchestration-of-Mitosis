Fluor_cell_combined = [Fluor_cell_01;Fluor_cell_02];
 Fluor_cell = Fluor_cell_combined;
pixel_size = 0.130;
per_pixels = 0.15;
time_pts_cut = 10;
time_pts_length = 5;
time_pts_cut_MCM = 15;
min_wt = 1.1;
max_wt = 2.0;

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
%figure (9)
num_cells = length(Fluor_cell(:,1));
peak_cyto = zeros(num_cells,4);
%% Protein 1(Syncut)

for i=1:num_cells
    
    
    cell_length_time = (cell2mat(Fluor_cell(i ,2)))*pixel_size;
    Cyto_Syn = cell2mat(Fluor_cell(i ,3));
    Nuc_Syn = cell2mat(Fluor_cell(i,4));
    
    CDK_activity_Syn_time = Nuc_Syn./ Cyto_Syn ;
    Cdc13_Cyto = cell2mat(Fluor_cell(i,5));
    Cdc13_Nuc = cell2mat(Fluor_cell(i,6));
    Total_cell = cell2mat(Fluor_cell(i, 7));
    CDK_Syn_smooth = smoothdata(CDK_activity_Syn_time,'sgolay',4);%default gaussian
    Cdc13_Cyto_smooth = smoothdata(Cdc13_Cyto,'sgolay',4);
    Cdc13_Nuc_smooth = smoothdata(Cdc13_Nuc, 'sgolay',4);
    Total_cell_smooth = smoothdata(Total_cell, 'sgolay',4);
%     [max_Cdc13, idx_Cdc13] = max(Cdc13_Nuc);
%     peak_cyto(i,1) = max(Cdc13_Cyto);
%     peak_cyto(i,2) = max(CDK_activity_Syn_time);
%     peak_cyto(i,3) = (max_Cdc13 - Cdc13_Nuc(idx_Cdc13 + 2));
%     peak_cyto (i,4) = max_Cdc13;
    time_vals = ((1:length(cell_length_time ))-1)*time_int;
    norm_time = time_vals/max(time_vals);
    
    figure(1)
    hold on

    plot(time_vals, CDK_activity_Syn_time);

    ylabel('N/C Ratio')
    xlabel('Time(minutes)')
    ylim([min_wt, max_wt])
    
    figure(2)
    hold on
    %yyaxis right
    %plot(time_vals, cell_length_time);
    %ylabel('Cell Length (um)')
    %yyaxis left
    plot(time_vals, Cdc13_Cyto);
   
    ylabel(' Mean Cdc 13 Cytoplasmic Intensity (a.u)')
    xlabel ('Time(minutes)')
    norm_Cdc13_Cyto = (Cdc13_Cyto-min(Cdc13_Cyto))/(max(Cdc13_Cyto)-min(Cdc13_Cyto));
    norm_Syn = (CDK_activity_Syn_time - min(CDK_activity_Syn_time))/(max(CDK_activity_Syn_time)-min(CDK_activity_Syn_time));

    norm_Cdc13_Cyto_smooth = (Cdc13_Cyto_smooth-min(Cdc13_Cyto_smooth))/(max(Cdc13_Cyto_smooth)-min(Cdc13_Cyto_smooth));
    norm_Syncut_smooth=(CDK_Syn_smooth - min(CDK_Syn_smooth))/(max(CDK_Syn_smooth)-min(CDK_Syn_smooth));
    figure(3)
    hold on
    subplot(1,2,2)
    
    plot(norm_Cdc13_Cyto_smooth, norm_Syncut_smooth)
    
    xlabel('Normalized Mean Cdc13 Cytoplasmic Intensity (a.u.)')
    ylabel('N/C Ratio Normalized')
    title(' Smoothed')
    hold on
    subplot(1,2,1)
    %plot(nuc_fluor_time_2/max(nuc_fluor_time_2), CDK_activity,'DisplayName',txt)
    plot(Cdc13_Cyto, CDK_activity_Syn_time)
    xlabel('Mean Cytoplasmic Cdc13 Intensity (a.u.)')
    ylabel('N/C Ratio')
    title('Raw Trace')
    figure (4)
    hold on
    plot(norm_Cdc13_Cyto_smooth,norm_Syncut_smooth)
    figure (5)
    hold on 
    plot(norm_time, CDK_activity_Syn_time)
    figure (6)
    hold on
    plot(norm_time, Cdc13_Nuc)
    figure (7)
    hold on
    plot(norm_time, Cdc13_Cyto);
    figure (8)
    hold on
    plot(norm_time, Total_cell)
    figure (9)
    hold on
    plot(Cdc13_Cyto_smooth,norm_Syncut_smooth)
    %text(max(time_vals_2),(nuc_fluor_time_2(end)),num2str(i))
    
    %text(max(time_vals_2),(nuc_fluor_time_2(end)),num2str(i))
end
hold off
legend show
%% 
figure(10)
num_cells_combined = 25;
for i=1:5:num_cells_combined
    cell_length_time = (cell2mat(Fluor_cell(i ,2)))*pixel_size;
    Cyto_Syn = cell2mat(Fluor_cell(i ,3));
    Nuc_Syn = cell2mat(Fluor_cell(i,4));
    
    CDK_activity_Syn_time = Nuc_Syn./ Cyto_Syn ;
    Cdc13_Cyto = cell2mat(Fluor_cell(i,5));
    Cdc13_Nuc = cell2mat(Fluor_cell(i,6));
    Total_cell = cell2mat(Fluor_cell(i, 7));
    CDK_Syn_smooth = smoothdata(CDK_activity_Syn_time,'sgolay',4);%default gaussian
    Cdc13_Cyto_smooth = smoothdata(Cdc13_Cyto,'sgolay',4);
    Cdc13_Nuc_smooth = smoothdata(Cdc13_Nuc, 'sgolay',4);
    Total_cell_smooth = smoothdata(Total_cell, 'sgolay',4);
    norm_Cdc13_Cyto_smooth = (Cdc13_Cyto_smooth-min(Cdc13_Cyto_smooth))/(max(Cdc13_Cyto_smooth)-min(Cdc13_Cyto_smooth));
    norm_Syncut_smooth=(CDK_Syn_smooth - min(CDK_Syn_smooth))/(max(CDK_Syn_smooth)-min(CDK_Syn_smooth));
    norm_syncut_raw = (CDK_activity_Syn_time - min(CDK_activity_Syn_time))/(max(CDK_activity_Syn_time)-min(CDK_activity_Syn_time));
    time_vals = ((1:length(cell_length_time ))-1)*time_int;
    norm_time = time_vals/max(time_vals);
    hold on
    plot(Cdc13_Cyto,norm_Syncut_smooth)
end
% elem_nonzeros = find(MCM_SYN(:,2));
% MCM_SYN = MCM_SYN(elem_nonzeros,:);
% Fluor_cell = Fluor_cell(elem_nonzeros,:);
% div_length = div_length(elem_nonzeros,:);
% div_time = div_time(elem_nonzeros,:);
% cov_Cdc13 = std(MCM_SYN(:,4))/mean(MCM_SYN(:,4));
% cov_SYN = std(MCM_SYN(:,5))/mean(MCM_SYN(:,5));
% MCM_SYN(1,16) = cov_Cdc13;
% MCM_SYN(1,17) = cov_SYN;