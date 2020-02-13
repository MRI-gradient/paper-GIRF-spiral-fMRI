% ROC plot analysis to evaluate spatial specificity of fMRI activation, as used in paper ?Feasibility of spiral fMRI based on an LTI gradient model?
% Nadine N Graedel, January 2020

% load data (masks and zstats)
maskGM = read_avw('mask_GM_bin_thr1p8_bin.nii.gz'); % GM mask = true positive mask
maskWM = read_avw('mask_WM_bin_thr1p8_bin.nii.gz'); % WM mask = false positive mask

zstat4moni=read_avw('imTot_moni_sli1to36_echo001_asl000_dyn1to100_dif000_swapdim_mxmypz.feat/stats/zstat4.nii.gz'); 
zstat4nomi=read_avw('imTot_nomi_k0corr_delaycorr_sli1to36_echo001_asl000_dyn1to100_dif000_swapdim_mxmypz.feat/stats/zstat4.nii.gz'); 
zstat4girf=read_avw('imTot_girf_k0corr_sli1to36_echo001_asl000_dyn1to100_dif000_swapdim_mxmypz.feat/stats/zstat4.nii.gz'); 

% determine number of true/false positives
dims = size(zstat4girf); 

zgirf = zstat4girf.*maskGM;
zmoni = zstat4moni.*maskGM;
znomi = zstat4nomi.*maskGM;

zgirfWM = zstat4girf.*maskWM;
zmoniWM = zstat4moni.*maskWM;
znomiWM = zstat4nomi.*maskWM;

zval = 2.3; 
zgirfact = zgirf(abs(zgirf) > zval); 
zmoniact = zmoni(abs(zmoni) > zval); 
znomiact = znomi(abs(znomi) > zval); 

zgirfactWM = zgirf(abs(zgirfWM) > zval); 
zmoniactWM = zmoni(abs(zmoniWM) > zval); 
znomiactWM = znomi(abs(znomiWM) > zval); 

disp(strcat('Number of true positives for GIRF = ', num2str(numel(zgirfact)))); 
disp(strcat('Number of true positives for Moni = ', num2str(numel(zmoniact)))); 
disp(strcat('Number of true positives for Nomi = ', num2str(numel(znomiact)))); 

disp(strcat('Number of false positives for GIRF = ', num2str(numel(zgirfactWM)))); 
disp(strcat('Number of false positives for Moni = ', num2str(numel(zmoniactWM)))); 
disp(strcat('Number of false positives for Nomi = ', num2str(numel(znomiactWM)))); 

disp(strcat('Average active z-stat for GIRF = ', num2str(mean(abs(zgirfact))))); 
disp(strcat('Average active z-stat for Moni = ', num2str(mean(abs(zmoniact))))); 
disp(strcat('Average active z-stat for Nomi = ', num2str(mean(abs(znomiact))))); 

disp(strcat('Maximum z-stat for GIRF = ', num2str(max(abs(zgirfact))))); 
disp(strcat('Maximum z-stat for Moni = ', num2str(max(abs(zmoniact))))); 
disp(strcat('Maximum z-stat for Nomi = ', num2str(max(abs(znomiact))))); 

disp(strcat('90th percentile z-stat for GIRF = ', num2str(prctile(abs(zgirfact),90)))); 
disp(strcat('90th percentile z-stat for Moni = ', num2str(prctile(abs(zmoniact),90)))); 
disp(strcat('90th percentile z-stat for Nomi = ', num2str(prctile(abs(znomiact),90)))); 

% Create ROC plot
zs = (0:0.1:10); 

for i = 1:numel(zs)
    zval = zs(i); 

    GMgirf(i) = numel(zgirf(abs(zgirf) > zval))/ numel(zgirf(abs(zgirf) > 0))*100; 
    GMmoni(i) = numel(zmoni(abs(zmoni) > zval))/numel(zmoni(abs(zmoni) > 0))*100; 
    GMnomi(i) = numel(znomi(abs(znomi) > zval))/numel(znomi(abs(znomi) > 0))*100; 

    WMgirf(i) = numel(zgirf(abs(zgirfWM) > zval))/numel(zgirf(abs(zgirfWM) > 0))*100; 
    WMmoni(i) = numel(zmoni(abs(zmoniWM) > zval))/numel(zmoni(abs(zmoniWM) > 0))*100; 
    WMnomi(i) = numel(znomi(abs(znomiWM) > zval))/numel(znomi(abs(znomiWM) > 0))*100; 

end 

figure; hold all;
grid on; 
plot(WMmoni,GMmoni,'-g','Linewidth',1.5); 
plot(WMgirf,GMgirf,'-b','Linewidth',1.5); 
plot(WMnomi,GMnomi,'-r','Linewidth',1.5); 
plot(0:0.1:100,0:0.1:100,'--k','Linewidth',1.5); 
lgd = legend('monitored','girf-predicted','nominal','Location','northwest'); 
lgd.FontSize = 14; 

zval = 2.3; 
GMgirf2p3 = numel(zgirf(abs(zgirf) > zval))/ numel(zgirf(abs(zgirf) > 0))*100; 
GMmoni2p3 = numel(zmoni(abs(zmoni) > zval))/numel(zmoni(abs(zmoni) > 0))*100; 
GMnomi2p3 = numel(znomi(abs(znomi) > zval))/numel(znomi(abs(znomi) > 0))*100; 

WMgirf2p3 = numel(zgirf(abs(zgirfWM) > zval))/numel(zgirf(abs(zgirfWM) > 0))*100; 
WMmoni2p3 = numel(zmoni(abs(zmoniWM) > zval))/numel(zmoni(abs(zmoniWM) > 0))*100; 
WMnomi2p3 = numel(znomi(abs(znomiWM) > zval))/numel(znomi(abs(znomiWM) > 0))*100;    

plot(WMgirf2p3,GMgirf2p3,'*b','Markersize',10, 'Linewidth',2); 
plot(WMmoni2p3,GMmoni2p3,'og','Markersize',10, 'Linewidth',2); 
plot(WMnomi2p3,GMnomi2p3,'+r','Markersize',10, 'Linewidth',2); 
xlabel('% of WM ROI voxels active','FontSize',14); 
ylabel('% of GM ROI voxels active','FontSize',14); 

% Output area under the curve
AUC_girf = trapz(fliplr(WMgirf),fliplr(GMgirf)) -100*100/2;
AUC_moni = trapz(fliplr(WMmoni),fliplr(GMmoni)) -100*100/2;
AUC_nomi = trapz(fliplr(WMnomi),fliplr(GMnomi)) -100*100/2;

disp(strcat('Area under the curve (AUC) for GIRF = ', num2str(AUC_girf))); 
disp(strcat('Area under the curve (AUC) for Moni = ', num2str(AUC_moni))); 
disp(strcat('Area under the curve (AUC) for Nomi = ', num2str(AUC_nomi))); 
