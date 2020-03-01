
%% GENERAL CASES

% Species with too many Phi<1
I_nanSub1=(Nobs_fracbelow1'>0.05);
% Species with too few data
I_nanNmin=(num_woa_unq'<5);

% For species with median habitat depth < zseas, use histograms from obs months
zseas=1e4;
I_z100=(Zhist(:,find(pct==90))<zseas & Zhist(:,find(pct==90))>0 & Ieo);
Phi_min(find(I_z100==1),:)=Phi_minM(find(I_z100==1),:);

%% SPECIAL CASES

% QC Phi_crit
Sp_cases{1}='Melanostigma pammelas'; % OMZ migrator with many Phi<1
Sp_cases{2}='Sergia tenuiremis'; % OMZ migrator with many Phi<1
Sp_cases{3}='Alitta succinea'; % benthic worm, poorly resolved habitat (esp. O2)
Sp_cases{4}='Bythograea thermydron'; % Hydrothermal vent crab, poorly resolved habitat (esp. T)
Sp_cases{5}='Carassius carassius'; % Crucian carp, freshwater fish
Sp_cases{6}='Oreochromis niloticus'; % Nile fish, largely freshwater
Sp_cases{7}='Pleuroncodes planipes'; % Tuna crab, targeted for O2 tolerance w/ v. few obs
Sp_cases{8}='Nematobrachion flexipes'; % Nematobrachion flexipes targeted for O2 tolerance
I_nanSp=Ieo*0;
for ii=1:length(Sp_cases)
    I_nanSp(find(strcmp(Sp_cases{ii},Par.Species_unq)))=1; 
end
I_nanSp(144)=1;

% QC for sensitivity to seasonal obs
clear Sp_cases
Sp_cases{1}='Acanthephyra purpurea'; % n=609
Sp_cases{2}='Acanthephyra smithi'; % n=7
Sp_cases{3}='Acipenser brevirostrum'; % n=29
Sp_cases{4}='Dosidicus gigas'; % 14
Sp_cases{5}='Gennadas valens'; % 
Sp_cases{6}='Penaeus vannamei'; % 
Sp_cases{7}='Sergia fulgens'; % 
I_nanFm=Phic_Fsspm==0;
for ii=1:length(Sp_cases)
    I_nanFm(find(strcmp(Sp_cases{ii},Par.Species_unq)))=1; 
end
Phic_Fsspm(find(I_nanFm==1))=nan;

Phic_Fssp(find(strcmp(Par.Species_unq,'Octopus vulgaris')))=nan;



%% set to NaN

I_nanPc = I_nanNmin | I_nanSp; 
Phi_min(find(I_nanPc==1),:)=nan;
Phi_minM(find(I_nanPc==1),:)=nan;

I_nanPc = I_nanFm | I_nanNmin | I_nanSp; 
Phic_F4d(find(I_nanPc==1))=nan;
Phic_Fssp(find(I_nanPc==1))=nan;
Phic_Fsspm(find(I_nanPc==1))=nan;

F_PhiT_4d(find(I_nanPc==1))=nan;
F_PhiO_4d(find(I_nanPc==1))=nan;
Fmax_ssp(find(I_nanPc==1))=nan;
Fmaxt_ssp(find(I_nanPc==1))=nan;
Fmaxo_ssp(find(I_nanPc==1))=nan;
Fmax_sspm(find(I_nanPc==1 | I_nanFm==1))=nan;
Fmaxt_sspm(find(I_nanPc==1 | I_nanFm==1))=nan;
Fmaxo_sspm(find(I_nanPc==1 | I_nanFm==1))=nan;

%%

% Species, #obs
Pobis.Species_unq=Par.Species_unq;
Pobis.num_tot_obis=num_tot_obis;
Pobis.num_woa_cell=num_woa_cell;
Pobis.num_woa_unq=num_woa_unq;
% habitat property histograms
Pobis.Phi_histWOA=Phi_histWOA;
Pobis.phi_wmax=nanmedian(Phi_max,2);
Pobis.phi_wmaxM=nanmedian(Phi_maxM,2);
Pobis.phi_wmin=nanmedian(Phi_min,2);
Pobis.phi_wminM=nanmedian(Phi_minM,2);
Pobis.Zhist=Zhist;
Pobis.Ohist=Ohist;
Pobis.Thist=Thist;
Pobis.T_atPatm=T_atPatm;
% F1 scores
Pobis.FratioT_4d=F_PhiT_4d;
Pobis.FratioO_4d=F_PhiO_4d;
Pobis.FratioT_ssp=Fmax_ssp./Fmaxt_ssp;
Pobis.FratioT_sspm=Fmax_sspm./Fmaxt_sspm;
Pobis.FratioO_ssp=Fmax_ssp./Fmaxo_ssp;
Pobis.FratioO_sspm=Fmax_sspm./Fmaxo_sspm;
% phi crit estimates
Pobis.phic_hist=nanmedian(Phi_min,2);
Pobis.phic_f4d=Phic_F4d';
Pobis.phic_fssp=Phic_Fssp';
Pobis.phic_fsspm=Phic_Fsspm';
%Pobis.phic_stsp=Phic_stsp(:,2);


%Pobis.phi_crit(strcmp('Cancer irroratus',Par.Species_unq),:) = 2.8;
