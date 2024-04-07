delimiterIn = ' ';
num_timepoints = 20;
time_intervals = 10;
multiplicity_thresh = 1;
phospho_dat = readtable('CCC6061_CCC6254_CCC6758 [21vv_imputation][5nn_smoothing] Phospho (STY)Sites_annotatedJC24 (2).csv');
localizations = importdata('JC_Localizations.mat');
%multi_array = table2array(phospho_dat(:,23));
% mult_find = find(contains(multi_array ,'1'));
% phospho_dat = phospho_dat(mult_find,:);
arr_phospho_all = 1./table2array (phospho_dat(:,1:20));

%arr_phospho = table2array(phospho_dat(:,28));
arr_sens = table2array(phospho_dat(:,40));
arr_gene_names = table2array(phospho_dat(:,38));
arr_prot_names = table2array(phospho_dat(:,36));
arr_prot_names = string(arr_prot_names);
joes_localization_names = table2array(localizations(1:end,1));
joes_localization_loc = string(table2array(localizations(1:end,2)));
joes_localization_names = string(joes_localization_names);
%arr_prot_names = num2str(cell2mat(arr_prot_names));
joe_filter = contains (arr_prot_names, joes_localization_names);


%%

arr_sens = arr_sens(joe_filter);
arr_phospho_all = arr_phospho_all(joe_filter,:);
%arr_phospho = arr_phospho(joe_filter,:);
arr_gene_names = arr_gene_names(joe_filter,:);
arr_prot_names = arr_prot_names(joe_filter,:);
num_sites = length (arr_sens(:,1));
%% 
arr_phosph_2 = strings(num_sites,1);
for site = 1:num_sites
    protein_select = arr_prot_names(site);
    joe_select = find(joes_localization_names == protein_select);
    arr_phosph_2(site) = joes_localization_loc(joe_select);

end
%% 

 idc_late = find(contains(arr_sens,'Late'));
 idc = strcmp(arr_phosph_2,'Cytoplasm'); 
 idc_SPB = strcmp(arr_phosph_2, 'SPB');
 idc_nuc = strcmp(arr_phosph_2, 'Nucleus');
 idc_NE = strcmp(arr_phosph_2, 'NE');
 idc_cyto_periphery = strcmp(arr_phosph_2, 'Cytoplasm (Periphery)');

%  idx_nuc = ~cellfun('isempty',idc_nuc);
%  idx_SPB = ~cellfun('isempty',idc_SPB);
%  idx = ~cellfun('isempty',idc);
cyto_elements = find(idc);
%cyto_elements(idx_SPB) = 0;
nuc_elements = find(idc_nuc);
SPB_elements = find (idc_SPB);
NE_elements = find (idc_NE);
cyto_periphery_elements = find(idc_cyto_periphery);

cyto_elements_2 = ismember(cyto_elements, idc_late);
cyto_elements = cyto_elements (cyto_elements_2);
nuc_elements_2 = ismember(nuc_elements,idc_late);
nuc_elements = nuc_elements (nuc_elements_2);
NE_elements_2 = ismember(NE_elements,idc_late);
NE_elements = NE_elements(NE_elements_2);
SPB_elements_2 = ismember(SPB_elements, idc_late);
SPB_elements = SPB_elements (SPB_elements_2);
cyto_periphery_elements_2 = ismember(cyto_periphery_elements,idc_late);
cyto_periphery_elements = cyto_periphery_elements(cyto_periphery_elements_2);

cyto_phosph = arr_phospho_all(cyto_elements,1:20); 
nuc_phosph = arr_phospho_all (nuc_elements,1:20);
SPB_phosph = arr_phospho_all (SPB_elements, 1:20);
NE_phosph = arr_phospho_all(NE_elements,1:20);
cyto_periphery_phosph = arr_phospho_all(cyto_periphery_elements,1:20);
arr_gene_cyt = arr_gene_names(cyto_elements);
arr_gene_nuc = arr_gene_names (nuc_elements);
arr_gene_SPB = arr_gene_names (SPB_elements);
arr_gene_NE = arr_gene_names (NE_elements);
arr_gene_cyt_peri = arr_gene_names(cyto_periphery_elements);
%% 

time_vals = [0,8,12,16,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130];

num_nuc_subs = length(nuc_phosph(:,1));


figure, 
for i = 1:num_nuc_subs

    plot(time_vals,nuc_phosph(i,:))
    text(time_vals(end), nuc_phosph(i,end), num2str(i))
    hold on
end
ylabel('Nuclear Phosphorylation')
xlabel('Time (Minutes)')
hold off
figure, 
num_cyto_subs = length(cyto_phosph(:,1));
for i = 1:num_cyto_subs
    plot(time_vals,cyto_phosph(i,:))
    text(time_vals(end), cyto_phosph(i,end), num2str(i))
    hold on
end
ylabel('Cytoplamsmic Phosphorylation')
xlabel('Time (Minutes)')
hold off

figure, 
num_NE_subs = length(NE_phosph(:,1));
for i = 1:num_NE_subs
    plot(time_vals,NE_phosph(i,:))
    text(time_vals(end), NE_phosph(i,end), arr_gene_NE(i))
    hold on
end
ylabel ('NE Phosphorylation')
xlabel ('Time (minutes)')
hold off

figure, 
num_cyto_per = length(cyto_periphery_phosph(:,1));
for i = 1:num_cyto_per 
    plot(time_vals,cyto_periphery_phosph(i,:))
    text(time_vals(end), cyto_periphery_phosph(i,end), num2str(i))
    hold on
end
ylabel ('Cytoplasmic Periphery Phosphorylation')
xlabel ('Time (minutes)')
hold off

figure, 
num_SPB = length(SPB_phosph(:,1));
for i = 1:num_SPB 
    plot(time_vals,SPB_phosph(i,:))
    text(time_vals(end), SPB_phosph(i,end), num2str(i))
    hold on
end
ylabel ('SPB Phosphorylation')
xlabel ('Time (minutes)')
hold off
% figure, 
% plot(arr_phospho_all(cyto_elements(6),1:20),'Color',[0 0.4470 0.7410],'LineStyle','-','LineWidth',1.0)
% hold on
% plot(arr_phospho_all(nuc_elements(4),1:20),'Color',[0.4660 0.6740 0.1880],'LineStyle','-','LineWidth',1.0)
% hold off

mean_cyto_phosph = mean(cyto_phosph);
mean_nuc_phosph = mean(nuc_phosph);
mean_SPB_phosph = mean(SPB_phosph);
mean_NE_phosph = mean(NE_phosph);
mean_cyto_periphery = mean(cyto_periphery_phosph);
med_cyto_phosph = median(cyto_phosph);
med_nuc_phosph = median(nuc_phosph);
med_SPB_phosph = median(SPB_phosph);
med_NE_phosph = median(NE_phosph);
med_cyto_periphery = mean(cyto_periphery_phosph);
 figure, 
 plot(time_vals,mean_nuc_phosph,'Color',[0.4660 0.6740 0.1880],'LineStyle','-','LineWidth',1.0) %green
 hold on
plot(time_vals, mean_cyto_phosph,'Color',[0 0.4470 0.7410],'LineStyle','-','LineWidth',1.0) %blue
plot(time_vals,mean_NE_phosph,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0)
plot(time_vals, mean_cyto_periphery,'Color',[0.6350 0.0780 0.1840],'LineStyle','-','LineWidth',1.0); 
plot (time_vals, mean_SPB_phosph, 'Color', [0.9290 0.6940 0.1250],'LineStyle','-','LineWidth',1.0);
ylabel('Mean Phosphorylation')
xlabel('Time (minutes)')
hold off



figure, 
 plot(time_vals,med_nuc_phosph,'Color',[0.4660 0.6740 0.1880],'LineStyle','-','LineWidth',1.0) %green
 hold on
plot(time_vals,med_cyto_phosph,'Color',[0 0.4470 0.7410],'LineStyle','-','LineWidth',1.0) %blue
plot(time_vals,med_NE_phosph, 'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',1.0)
plot(time_vals, med_cyto_periphery,'Color',[0.6350 0.0780 0.1840],'LineStyle','-','LineWidth',1.0); 
plot (time_vals, med_SPB_phosph, 'Color', [0.9290 0.6940 0.1250],'LineStyle','-','LineWidth',1.0);
ylabel('Median Phosphorylation')
xlabel('Time (minutes)')
hold off
%% 

all_phosph_late = arr_phospho_all(idc,:);
figure, 
plot(time_vals, mean(all_phosph_late))
