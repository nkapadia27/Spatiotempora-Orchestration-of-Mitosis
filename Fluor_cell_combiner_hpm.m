%Fluor_cell_combined = [Fluor_cell_01;Fluor_cell_02];
Fluor_cell = Fluor_cell_combined;
pixel_size = 0.130;
per_pixels = 0.15;
time_pts_cut = 10;
time_pts_length = 5;
time_pts_cut_MCM = 15;
time_int = 10; %minutes
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
    
   cell_length_time = (cell2mat(Fluor_cell(i ,2)))*pixel_size;
   Cyto_MCM = cell2mat(Fluor_cell(i,3));
   Nuc_MCM = cell2mat(Fluor_cell(i,4));
    
    CDK_activity_MCM_time = Cyto_MCM./ Nuc_MCM ;
    Cdc13_Cyto = cell2mat(Fluor_cell(i,5));
    Cdc13_Nuc = cell2mat(Fluor_cell(i,6));
    Total_cell = cell2mat(Fluor_cell(i, 7));
    CDK_MCM_smooth = smoothdata(CDK_activity_MCM_time,'sgolay',5);
%     [pks,locs] = findpeaks(CDK_MCM_smooth,'Threshold',0.005);
%     id_MCM = locs;
    Cdc13_Cyto_smooth = smoothdata(Cdc13_Cyto,'sgolay',4);
    Cdc13_Nuc_smooth = smoothdata(Cdc13_Nuc,'sgolay',4);
    Total_cell_smooth = smoothdata(Total_cell, 'sgolay',4);
    norm_Cdc13_Nuc_raw = (Cdc13_Nuc-min(Cdc13_Nuc))/(max(Cdc13_Nuc)-min(Cdc13_Nuc));
    norm_MCM_raw = (CDK_activity_MCM_time-min(CDK_activity_MCM_time))/(max(CDK_activity_MCM_time)-min(CDK_activity_MCM_time));
    
    time_vals = ((1:length(cell_length_time ))-1)*time_int;
    %time_vals_recentered = ((1-id_MCM):length(cell_length_time)-id_MCM)*time_int;
    norm_time = time_vals/max(time_vals);
    
    figure(1)
    hold on
    plot(time_vals, CDK_activity_MCM_time);
    ylabel('MCM sensor readout')
    xlabel('Time(minutes)')

    figure(2)
    hold on

    plot(time_vals, Cdc13_Nuc);
    text(max(time_vals),(Cdc13_Nuc(end)),num2str(i))
    ylabel(' Mean Cdc 13 Nuclear Intensity (a.u)')
    xlabel ('Time(minutes)')
    
    figure(3)
    hold on
    subplot(1,2,2)
    norm_Cdc13_nuc = (Cdc13_Nuc_smooth-min(Cdc13_Nuc_smooth))/(max(Cdc13_Nuc_smooth)-min(Cdc13_Nuc_smooth));
    norm_MCM = (CDK_MCM_smooth-min(CDK_MCM_smooth))/(max(CDK_MCM_smooth)-min(CDK_MCM_smooth));
    plot(norm_Cdc13_nuc,norm_MCM)
    xlabel('Mean Normalized Nuclear Cdc13 Intensity (a.u.)')
    ylabel('C/N Ratio Normalized')
    title('Normalized Smoothed Trace')
    hold on
    subplot(1,2,1)
    plot(Cdc13_Nuc, CDK_activity_MCM_time)
    xlabel('Mean Nuclear Cdc13 Intensity (a.u.)')
    ylabel('C/N Ratio')
    title('Raw Trace')
    figure (4)
    hold on
    plot(norm_Cdc13_nuc,norm_MCM)
    figure (5)
    hold on 
    plot(time_vals, CDK_MCM_smooth)
    figure (6)
    hold on
    plot(time_vals, Cdc13_Nuc)
    figure (7)
    hold on
    plot(norm_time, Cdc13_Cyto);
    figure (8)
    hold on
    plot(norm_time, Total_cell)
    figure (9)
    hold on

    plot(Cdc13_Nuc_smooth, norm_MCM)
    %text(max(time_vals_2),(nuc_fluor_time_2(end)),num2str(i))
end
hold off

% elem_nonzeros = find(MCM_SYN(:,2));
% MCM_SYN = MCM_SYN(elem_nonzeros,:);
% Fluor_cell = Fluor_cell(elem_nonzeros,:);
% div_length = div_length(elem_nonzeros,:);
% div_time = div_time(elem_nonzeros,:);
% cov_Cdc13 = std(MCM_SYN(:,4))/mean(MCM_SYN(:,4));
% cov_MCM = std(MCM_SYN(:,5))/mean(MCM_SYN(:,5));
% MCM_SYN(1,16) = cov_Cdc13;
% MCM_SYN(1,17) = cov_MCM;
% num_cells = length(Fluor_cell(:,1));