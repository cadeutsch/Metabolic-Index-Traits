
%% GENERAL CASES

% Species with too many Phi<1
I_nanSub1=(Nobs_fracbelow1'>0.05);
% Species with too few data
I_nanNmin=(num_woa_unq'<5);


%% SPECIAL CASES

% QC Phi_crit
% 'Alitta succinea' benthic worm, poorly resolved habitat (esp. O2)
% 'Bythograea thermydron' Hydrothermal vent crab, poorly resolved habitat (esp. T)
% 'Carassius carassius' Crucian carp, freshwater fish
% 'Oreochromis niloticus' Nile fish, largely freshwater
% 'Pleuroncodes planipes' Tuna crab, targeted for O2 tolerance, poorly sampled
% 'Nematobrachion flexipes' targeted for O2 tolerance
Sp_cases={'Alitta succinea','Bythograea thermydron','Carassius carassius',...
    'Oreochromis niloticus','Pleuroncodes planipes','Nematobrachion flexipes'};
I_nanSp=Ieo*0;
for ii=1:length(Sp_cases)
    I_nanSp(find(strcmp(Sp_cases{ii},Par.Species_unq)))=1; 
end

% set to NaN
I_nanPc = I_nanNmin | I_nanSp; 
Phi_min(find(I_nanPc==1),:)=nan;
Phi_minM(find(I_nanPc==1),:)=nan;
Phic_F4d(find(I_nanPc==1),:)=nan;
Phic_Fssp(find(I_nanPc==1),:)=nan;
Phic_Fsspm(find(I_nanPc==1),:)=nan;
F_PhiT_4d(find(I_nanPc==1),:)=nan;
F_PhiO_4d(find(I_nanPc==1),:)=nan;
Fmax_ssp(find(I_nanPc==1),:)=nan;
Fmaxt_ssp(find(I_nanPc==1),:)=nan;
Fmaxo_ssp(find(I_nanPc==1),:)=nan;
Fmax_sspm(find(I_nanPc==1),:)=nan;
Fmaxt_sspm(find(I_nanPc==1),:)=nan;
Fmaxo_sspm(find(I_nanPc==1),:)=nan;

%%

% use all months when <30 grid points
I_lowN=find(num_woa_unq<30 & num_woa_unq>0);
Phic_fssp_merge=Phic_Fsspm;
Phic_fssp_merge(I_lowN,:)=Phic_Fssp(I_lowN,:);
Phi_min_merge=Phi_minM;
Phi_min_merge(I_lowN,:)=Phi_min(I_lowN,:);
Phi_min0_merge=Phi_minM0;
Phi_min0_merge(I_lowN,:)=Phi_min0(I_lowN,:);
Phi_max_merge=Phi_maxM;
Phi_max_merge(I_lowN,:)=Phi_max(I_lowN,:);

% use 5% hist when 5/10 are >dmax apart
dmax=0.5;
d5=abs(diff(Phi_min_merge,[],2));
df=abs(Phic_F4d(:,1)-nanmedian(Phi_min_merge,2));
I_d5=find(d5>dmax & df>dmax);
Phi_min_merge(I_d5,2)=nan;
Phi_min0_merge(I_d5,2)=nan;

%%  Set phi_crit

% Baseline 
phi_crit=nanmedian(Phi_min_merge(:,1:2),2);

% Individual species adjustments
I=find(strcmp(Pobis.Species_unq,'Chromis atripectoralis')); 
phi_crit(I)=nanmedian(Phic_fssp_merge(I,1),2);
I=find(strcmp(Pobis.Species_unq,'Gadus morhua')); 
phi_crit(I)=nanmedian(Phic_fssp_merge(I,1),2);
Phic_F4d(I,1)=nanmedian(Phic_fssp_merge(I,1),2);
%Pobis.phic_f4d(I,1)=median([Phic_fssp_merge(I,2) Phic_fssp_mergem(I,2)]);
I=find(strcmp(Pobis.Species_unq,'Morone saxatilis'));
phi_crit(I)=nanmedian(Phic_fssp_merge(I,1),2);
I=find(strcmp(Pobis.Species_unq,'Notostomus gibbosus'));
phi_crit(I)=Phic_F4d(I,1);
I=find(strcmp(Pobis.Species_unq,'Sciaenops ocellatus'));
phi_crit(I)=nanmedian(Phic_fssp_merge(I,1),2);
phi_crit([141],:)=nan; % Atl. cod with LC50 (not Pcrit)
%phi_crit(strcmp('Cancer irroratus',Par.Species_unq),:) = 2.8;

%% Save to structure (Pobis)

% Species, #obs
Pobis.Species_unq=Par.Species_unq;
Pobis.num_tot_obis=num_tot_obis;
Pobis.num_woa_cell=num_woa_cell;
Pobis.num_woa_unq=num_woa_unq;
Pobis.Nobs_fracbelow1=Nobs_fracbelow1;
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
Pobis.phic_hist=Phi_min_merge;
Pobis.phic_f4d=Phic_F4d;
Pobis.phic_fssp=Phic_fssp_merge;
Pobis.phi_crit=phi_crit;

