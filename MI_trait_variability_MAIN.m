
clear all
close all

set(groot,'defaultLineLineWidth', 2)
set(groot,'DefaultAxesFontSize', 14)
set(groot,'DefaultLegendFontSize', 14)

%% Set constants

global kb Rgas T0 Tref Patm atm2kPa

kb=8.6e-5; 
Rgas=8.3e-3;
T0=273.15;
Tref=15;
Patm=0.209;
atm2kPa=101.325;

%% Load data 

% Laboratory Data
    
load Dat.mat

Species_unq=unique(Dat.Species);
Dat.Species_unq=Species_unq;
nSp=length(Species_unq);

% Metabolic Index parameters 

load Par.mat
    
% Metabolic rate params
del_norm=0.25; % intra-specific B normalization (0 or 0.25)
alphaD=Par.M.est(:,1).*Par.B.^del_norm; 
alphaD_altB=Par.M.est(:,1).*Par.B.^0.25; 
del=Par.M.est(:,2); 
Emet=Par.M.est(:,3); 
alphaD_err=Par.M.err(:,1); 
Emet_err=Par.M.err(:,3); 

% Hypoxia (Met.Idx.) params
Ao=Par.P.est(:,1); 
eps=Par.P.est(:,2); 
Eo=Par.P.est(:,3); 
Ao_err=Par.P.err(:,1); 
Eo_err=Par.P.err(:,3); 
Esup=Emet-Eo;
alphaS = alphaD.*Ao;
alphaSbar=nanmean(alphaS);
tmp=((alphaD_err./alphaD).^2 + (Ao_err./Ao).^2).^0.5;
alphaS_err = tmp.*alphaS;

% Flags & counts
Ieo=Par.Ieo;
Iem=Par.Iem;
In2=Par.In2;
nSp_Eo_all=sum(Ieo);
nSp_Eo_sig=length(find(Par.P.pval(:,3)<.05));
Isig=(find(Par.P.pval(:,3)<.05));

%% Additional Data

% Categorize Taxa
phynames={'Chordata','Mollusca','Crustacea','Cnidaria','Tunicata','Annelida','Other'};
Par = taxon_categ(Par,phynames);
Dat = taxon_categ(Dat,phynames);

% Ecological body mass
Par.B(isnan(Par.B))=1;
Par.Beco=Par.B;     
Par.Beco(find(strcmp('Zoarces viviparus',Par.Species_unq)))=300;
Par.Beco(find(strcmp('Diplodus puntazzo',Par.Species_unq)))=1000;
% C. atrip. based on adult mass of C. viridis scaled by their length ratio 12/10
Par.Beco(find(strcmp('Chromis atripectoralis',Par.Species_unq)))=5*(12/10)^3; 

% Thermal tolerance data
load Tol
Imar=Tol.Imar;
Itmc=Tol.Itmc;
Itax=Imar;
CTmax=ones(length(Par.B),1)*nan;
for i=1:length(Par.Species_unq)
    I=find(strcmp(Par.Species_unq(i),Tol.Sp_name));
    if ~isempty(I)
        CTmax(i)=Tol.Tmax(I);
    end
end

% Eo variation
load Eo_TempLoHi.mat
Icw=isfinite(ParC.Eo-ParW.Eo); Icw(46)=0; % second cod value
dEdTbar(1)=1e-4; % No Eo trend 
dEdTbar(2)=-nanmean(Par.P.param_EoT(:,4)); % Intraspecific mean #1
dEdTbar(3)=nanmean(ParW.Eo(Icw)-ParC.Eo(Icw))/15; % Intraspecific mean #2
dEdTbar(4)=(nanmean(ParW.Eo)-nanmean(ParC.Eo))/15; % Intraspecific mean #2

% Get Egas, Event, Ecirc
load Es
% SMS+FAS (Marine & Terrestrial)
load FAS

%% Global T/O2 Fields
        
load WOA

% Hydrographic grid
dT=1; Tedge=-4:dT:34; Tcent=Tedge(2:end)-dT/2;
dO=0.01; Oedge=0:dO:0.25; Ocent=Oedge(2:end)-dO/2;

% T/O volume histogram
vol4=repmat(WOA.v3d,[1 1 1 12]);
[bin] = bindata(WOA.temp(:),Tedge,WOA.po2(:),Oedge,vol4(:)./vol4(:),vol4(:));
TOvol=bin.P;

%% Biogeographic metrics (OBIS)

pct=[0 1 5:5:95 99 100]; 
pctmin=[5 10];
tmin=10;
brange=[0.5:.25:1];

try
%    uncomment this line to recompute Pobis
    load Pobis.mat
catch
    Phi_min=ones(length(Par.B),length(pctmin))*nan;
    Phi_max=ones(length(Par.B),length(pctmin))*nan;
    F_PhiT=ones(length(Par.B),1)*nan;F_PhiT_4d=F_PhiT;
    F_PhiO=ones(length(Par.B),1)*nan;F_PhiO_4d=F_PhiO;
    Nobs_fracbelow1=ones(length(Par.B),1)*nan;
    Thist=ones(length(Par.B),length(pct))*nan; Ohist=Thist; PO2hist=Thist; Zhist=Thist;
    Tmax=ones(length(Par.B),length(pct))*nan; T_atPatm=Tmax;
    Phic_Fssp=ones(length(Par.B),length(brange))*nan;Phic_Fsspm=Phic_Fssp; Phic_F4d=Phic_Fssp;
    Fmax_ssp=Phic_Fssp; Fmax_sspm=Phic_Fssp;
    Fmaxt_ssp=Phic_Fssp; Fmaxt_sspm=Phic_Fssp;
    Fmaxo_ssp=Phic_Fssp; Fmaxo_sspm=Phic_Fssp;
    for i=1:length(Par.Species_unq)
        spec1=Par.Species_unq(i);
        fname=strrep(spec1{1},' ','_')
        [Pdat] = extract_Spec_MIparm(Par,spec1);
        if isfinite(Pdat.Eo) && i~=141 % species has physiological data
            [Odat] = get_obis_data(spec1,WOA,'auto',obis_dir);
            if ~isempty(Odat.depth)  % species has distribution data
                num_tot_obis(i)=Odat.num_tot_obis;
                num_woa_cell(i)=Odat.num_woa_cell;
                num_woa_unq(i)=Odat.num_woa_unq;
                                    
                load(['output_OBIS/output_' fname '_mod'])

                % histogram analysis
                [C,IA,IB]=intersect(Pchis.pct,pctmin);
                Phi_min(i,:) = Pchis.Phi_decile(IA);
                Phi_minM(i,:) = Pchis.Phi_decileM(IA);
                [C,IA,IB]=intersect(Pchis.pct,100-pctmin);
                Phi_max(i,:) = Pchis.Phi_decile(IA);
                Phi_maxM(i,:) = Pchis.Phi_decileM(IA);
                Phi_histWOA(i,:) = Pchis.Phi_decileWOA;
                Thist(i,:)=Pchis.T_decile;
                Ohist(i,:)=Pchis.O2_decile;
                PO2hist(i,:)=Pchis.pO2_decile;
                Zhist(i,:)=prctile(Odat.depth(Odat.iz==1),pct);
                Nobs_fracbelow1(i)=sum(Pchis.PhistNobs(1:2))/sum(Pchis.PhistNobs);
                T_atPatm(i,:)=Pchis.T_atPatm';
                
                % 4-dimensional analysis
                [Phic_F4d(i), F_PhiT_4d(i), F_PhiO_4d(i)]=get_Fstats(Pc4d);
                
                % State-space analysis
                if ~isempty(Pcssp) && ~strcmp(spec1,'Bythograea thermydron')
                    [Phic_stsp(i,:)] = get_StSp( Pcssp , Pdat, 0);
                    [Fmax_ssp(i),idx]=max(Pcssp.Fssp.f);
                    Phic_Fssp(i)=Pcssp.Fssp.p(idx);
                    [Fmax_sspm(i),idx]=max(Pcssp.Fsspm.f);
                    Phic_Fsspm(i)=Pcssp.Fsspm.p(idx);
                    [Fmaxt_ssp(i),idx]=nanmax(Pcssp.Fsst.f(Pcssp.Fsst.p>Pchis.T_decile(12)));
                    [Fmaxt_sspm(i),idx]=nanmax(Pcssp.Fsstm.f(Pcssp.Fsstm.p>Pchis.T_decile(12)));
                    [Fmaxo_ssp(i),idx]=nanmax(Pcssp.Fsso.f(Pcssp.Fsso.p<Pchis.pO2_decile(12)));
                    [Fmaxo_sspm(i),idx]=nanmax(Pcssp.Fssom.f(Pcssp.Fssom.p<Pchis.pO2_decile(12)));
                    
                end
            end
        end
    end
    
    filter_eTraits
    
    save Pobis.mat Pobis
    
end


%% ATmax

clear Tc Ta Ta2
for i=1:length(dEdTbar)
Par.dEdT=ones(size(Par.B))*dEdTbar(i); 
[ Tc(:,i), Ta(:,i), Ta2(:,i) ] = compute_ATmax(Par,Pobis);
end
% Species-specific Eo(T) where available
Par.dEdT=-Par.P.param_EoT(:,4); 
I=isnan(Par.dEdT); Par.dEdT(I)=0.01;
[ Tc(:,i+1), Ta(:,i+1), Ta2(:,i+1) ] = compute_ATmax(Par,Pobis); 
Tc(I,i+1)=nan; Ta(I,i+1)=nan; Ta2(I,i+1)=nan;

% find species with both phi_crit and dEdT estimates
I=isfinite(Ta(:,5)); 
% replace their estimates based on mean-species parameters
for i=2:4; 
    Ta(I,i)=Ta(I,5); 
end
% find species without phi_crit estimates
I=isnan(Ta); 
% replace their estimates based on mean-species parameters
Ta_combo=Ta; Ta_combo(I)=Ta2(I); 
%Ta_combo(144,:)=nan;

%% Make plots

% Main figs
MI_traits_Figs




