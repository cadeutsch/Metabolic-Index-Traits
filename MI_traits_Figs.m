set(gcf,'defaultLineLineWidth', 2)
set(gcf,'DefaultAxesFontSize', 12)
set(gcf,'DefaultLegendFontSize', 14)

doPrint=0;

clear frq_*

%% Manuscript Fig1

set(gcf,'defaultLineLineWidth', 2)
set(gcf,'DefaultAxesFontSize', 13)
set(gcf,'DefaultLegendFontSize', 14)

figure(1); clf reset;
wysiwyg

dV=1; Vedges=0:dV:30; 
dM=0.3; Medges=-2:dM:4; 
dE=0.15; Eedges=-1:dE:2; 
Tedges=0:2:50; 

utaxM=phynames; 

for i=1:length(utaxM)
    Itax=strcmp(utaxM{i},Par.categ);
    frq_alphaD(:,i) = histc(log10(alphaD(Iem&Itax)),Medges);
    frq_Emet(:,i) = histc(Emet(Iem&Itax),Eedges);
end

utaxP=utaxM; 
for i=1:length(utaxP)
    frq_Vh(:,i) = histc((100./Ao(strcmp(utaxP{i},Par.categ))),Vedges);
    frq_Eo(:,i) = histc((Eo(strcmp(utaxP{i},Par.categ))),Eedges);
end

% Plotting

[cmap]=cbrewer('qual','Set1',length(utaxM));
colormap(cmap)
msize=50;
fsize=10;

figure(1); clf reset;

subplot('position',[0.25 0.76 0.25 0.22])
Tstar=1/kb*(1./(Tedges+T0) - 1/(Tref+T0));
pcrit=.04*exp(-0.4 * Tstar);
plot(Tedges,pcrit,'b-'); hold on;
plot(Tedges,pcrit*3,'r')
plot(Tedges,Tedges*0+Patm,'k-')
tlim=[25.2 47.8]; ylim([0 Patm*1.2])
plot([tlim(1) tlim(1)],[0 Patm],'r--')
plot([tlim(2) tlim(2)],[0 Patm],'b--')
plot([Tref Tref],[0 Patm*0.19],'b:')
plot([0 Tref],[Patm Patm]*0.19,'b:')
plot([Tref Tref],[0.19 0.57]*Patm,'r:')
plot([0 Tref],[Patm Patm]*0.57,'r:')
set(gca,'XTick',[], 'YTick', [])

% panel A
subplot('position',[0.13 0.55 0.24 0.17]) % [left bottom width height]
bar(Medges+dM/2,frq_alphaD,'stacked')
ylabel('# Species'); xlabel('log10 Resting Metabolic Rate (\alpha_D)');
axis([-3 4 0 45])
h2=legend(utaxM,'Location','NorthWest');
set(h2,'FontSize',8)

% panel B
subplot('position',[0.42 0.55 0.24 0.17]) % [left bottom width height]
bar(Eedges+dE/2,frq_Emet,'stacked')
xlabel('Temp. sensitivity, E_{d} [eV]'); ylabel('# Species'); 
axis([-.5 2.5 0 45])
h2=legend(utaxM,'Location','NorthEast');
set(h2,'FontSize',8)

% panel C
subplot('position',[0.13 0.32 0.24 0.17]) % [left bottom width height]
bar(1e-2*(Vedges+dV/2),frq_Vh,'stacked')
ylabel('# Species');xlabel('Hypoxia Vulnerability (V_h) [atm]');
axis([-0.01 0.2 0 21])
hold on;

% panel D
subplot('position',[0.42 0.32 0.24 0.17]) % [left bottom width height]
bar(Eedges+dE/2,frq_Eo,'stacked')
xlabel('Temp. sensitivity, E_{o} [eV]'); ylabel('# Species'); 
axis([-.5 2.5 0 21])

% panel E
subplot('position',[0.13 0.06 0.24 0.2]) % [left bottom width height]
gV=[0.02 0.05 0.1];
h1=plot([0 30],[0 30]/gV(1),'r--'); hold on
t1=text(1,450,['V_h = ' num2str(gV(1))],'Color','r','FontSize',12)
h2=plot([0 30],[0 30]/gV(2),'g--'); hold on
t2=text(15,270,['V_h = ' num2str(gV(2))],'Color','g','FontSize',12)
h3=plot([0 30],[0 30]/gV(3),'b--'); hold on
t3=text(20,150,['V_h = ' num2str(gV(3))],'Color','b','FontSize',12)
I=find(alphaS>0); 
Icol=Par.categ_num(I);
scatter(alphaD(I),alphaS(I),Icol./Icol*msize,cmap(Icol,:),'filled'); hold on;
H=errorbarxy(alphaD(I),alphaS(I), alphaD_err(I),alphaS_err(I), {'k', 'k', 'k'});
H.hMain.LineStyle='none';
x=alphaD(I); y=alphaS(I);
[B,BINT,R,RINT,STATS] = regress(y,[x x./x]);
xlim([0 30]); ylim([0 500]);
ylabel('\alpha_S');xlabel('\alpha_D')

% panel F
subplot('position',[0.42 0.06 0.24 0.2]) % [left bottom width height]
xp=[-2,2,2,-2]; 
transp=0.2;
y1=Es.bioCold(3); dy=Es.bioCold(4); yp=xp; yp(1:2)=yp(1:2)-(y1); yp(3:4)=yp(3:4)-(dy);
h1=fill(xp,yp,'blue','FaceAlpha',transp,'EdgeColor',[1 1 1]); hold on;
y1=Es.bioWarm(3); dy=Es.bioWarm(4); yp=xp; yp(1:2)=yp(1:2)-(y1); yp(3:4)=yp(3:4)-(dy);
h2=fill(xp,yp,'red','FaceAlpha',transp,'EdgeColor',[1 1 1]);
yp=xp; yp(1:2)=yp(1:2)-(0.42); yp(3:4)=yp(3:4)-(0.21);
h3=fill(xp,yp,'yellow','FaceAlpha',0.4,'EdgeColor','yellow','EdgeAlpha',0.4);

scatter(Emet(I),Eo(I),Icol./Icol*msize,cmap(Icol,:),'filled'); hold on;
plot(-1:2,-1:2,'k-'); 
x=Emet(I); y=Eo(I);
[B,BINT,R,RINT,STATS] = regress(y,[x x./x]);
plot([-.5 2],[-.5 2]*B(1)+B(2),'k:')
ylabel('E_o [eV]');xlabel('E_{d} [eV]'); 
H=errorbarxy(Emet(I),Eo(I), Par.M.err(I,3),Par.P.err(I,3), {'k', 'k', 'k'});
H.hMain.LineStyle='none';
xlim([-.5 2]); ylim([-.5 1.5]);

x=alphaD(I); y=alphaS(I); [B,BINT,R,RINT,STATS] = regress(y,[x x./x]);
x=alphaD(Isig); y=alphaS(Isig); [B,BINT,R,RINT,STATS] = regress(y,[x x./x]);

x=Emet(I); y=Eo(I); [B,BINT,R,RINT,STATS] = regress(y,[x x./x]);
x=Emet(Isig); y=Eo(Isig); [B,BINT,R,RINT,STATS] = regress(y,[x x./x]);

std(100./Ao(isfinite(Ao)))/mean(100./Ao(isfinite(Ao)))
iqr(100./Ao(isfinite(Ao)))/median(100./Ao(isfinite(Ao)))
std(alphaD(Iem&Itax))/mean(alphaD(Iem&Itax))
iqr(alphaD(Iem&Itax))/median(alphaD(Iem&Itax))

%if doPrint
print -djpeg MI_traits_Fig1.jpg
print(gcf,'-depsc2', '-painters', 'MI_traits_Fig1.eps');
%end

%% Manuscript Fig3

[cmap]=cbrewer('seq','YlGnBu',11);
cmap2=cmap(5:end,:);
cmap2(1,:)=[1 1 1]*0.5;
colormap(cmap2)

xl=[8 33];
cax=[-.5 3];

figure(3); clf reset;

subplot(231); i=find(strcmp(Par.Species_unq,'Paralichthys dentatus')); pc_F3=Pobis.phic_fssp(i); dopanel; 
[Par.P.est(i,3) Pobis.phic_hist(i) Pobis.phic_fssp(i) Pobis.phic_f4d(i)]
subplot(232); i=find(strcmp(Par.Species_unq,'Nautilus pompilius')); pc_F3=Pobis.phic_fssp(i); dopanel; 
[Par.P.est(i,3) Pobis.phic_hist(i) Pobis.phic_fssp(i) Pobis.phic_f4d(i)]
subplot(233); i=find(strcmp(Par.Species_unq,'Styela plicata')); pc_F3=Pobis.phic_fssp(i); dopanel; 
[Par.P.est(i,3) Pobis.phic_hist(i) Pobis.phic_fssp(i) Pobis.phic_f4d(i)]

if doPrint
print(gcf,'-djpeg','MI_traits_Fig3.jpg')
print(gcf,'-depsc2', '-painters', 'MI_traits_Fig3.eps');
end


%% Manuscript Fig4

clear frq_*

dE=1;
eoffset=dE/2;
w=0.5;ht=0.22;
loc='NorthEast'; %loc='NorthEast';
m2s=0.5;

ncut=Pobis.num_woa_cell;
nmin=50; mincat=['n < ' num2str(nmin)];
Inum=Pobis.num_woa_cell>nmin;
categ2=Par.categ;
for i=1:length(categ2)
if Inum(i)==0; categ2{i}=mincat; end
end

edges=0:dE:10; 
O2edges=0:5:400;
pO2edges=0:0.01:0.25;
utaxM=unique(categ2);
for i=1:length(utaxM)
    Itax=strcmp(utaxM{i},categ2);
    frq_phic(:,i) = histc(Pobis.phi_crit(Ieo&Itax),edges);
    frq_tmax(:,i) = histc(Pobis.Thist(Ieo&Itax,end-1),Tedge);
    frq_tmin(:,i) = histc(Pobis.Thist(Ieo&Itax,2),Tedge);
    frq_omin(:,i) = histc(Pobis.Ohist(Ieo&Itax,2),O2edges);
%    frq_po2min(:,i) = histc(Pobis.PO2hist(Ieo&Itax,2),pO2edges);
end

Fedges=0:20;
mms=mmsVal; mms(mms==0)=nan;
utaxA=unique(mmsTax);
for i=1:length(utaxA)
    Itax=strcmp(utaxA{i},mmsTax);
    frq_mms(:,i) = histc(m2s*(mms(Itax)+1),edges);
    frq_fas(:,i) = histc(mms(Itax),Fedges);
end
mmsTax=mmsTax(sum(frq_mms)>0);
frq_mms=frq_mms(:,sum(frq_mms)>0);
frq_fas=frq_fas(:,sum(frq_fas)>0);

nsT=size(SMS,1);
nbinT=ones(nsT,length(edges))*nan;
for i=1:nsT
nbinT(i,:) = histc(squeeze(SMS(i,:)), edges);
end

% Plot

fsize=12;
set(gcf,'DefaultAxesFontSize', 16)

% top
figure(4);clf reset;
[cmap]=cbrewer('qual','Set1',length(utaxM));
cmap=cat(1,cmap,[1 1 1]*0.95); % gray bar at the end
colormap(cmap)
subplot('Position',[0.25 0.7 w ht]); 
bar(eoffset+edges,frq_phic,'stacked');
h1=legend(utaxM,'Location',loc);
axis([0 10 0 25])
xlabel('Metabolic Habitat Limit (\Phi^{crit})'); ylabel('# Species'); 
set(h1,'FontSize',fsize)
if doPrint
print -djpeg MI_traits_Fig4a.jpg
print(gcf,'-depsc2', '-painters', 'MI_traits_Fig4a.eps');
end

% middle
figure(4);clf reset;
[cmap]=cbrewer('qual','Set1',9);
colormap(cmap)
subplot('Position',[0.25 0.4 w ht]); 
bar(eoffset+edges,frq_mms,'stacked')
h2=legend(utaxA,'Location','NorthEast'); set(h2,'FontSize',fsize-2)
xlabel('Sustained Metabolic Scope (SMS)');ylabel('# Species'); 
axis([0 10 0 35])
if doPrint
print -djpeg MI_traits_Fig4b.jpg
print(gcf,'-depsc2', '-painters', 'MI_traits_Fig4b.eps');
end

% bottom
figure(4);clf reset;
subplot('Position',[0.25 0.1 w ht]); 
bar(eoffset+edges,0.5*nbinT','stacked'); 
xlabel('Sustained Metabolic Scope (SMS)');ylabel('# Species'); 
h3=legend(ttaxa,'FontSize',10,'Location',loc);  
axis([0 10 0 20])
set(h3,'FontSize',fsize)
[cmap]=cbrewer('qual','Set1',4);
colormap(cmap)
if doPrint
print -djpeg MI_traits_Fig4c.jpg
print(gcf,'-depsc2', '-painters', 'MI_traits_Fig4c.eps');
end


%% Manuscript Fig5

Tedges=0:2:50; Tcent=Tedges(2:end)-mean(diff(Tedges))/2;
Tedges1=1:2:51; Tcent1=Tedges1(2:end)-mean(diff(Tedges1))/2;
Tcphi = histc(Tc(:,4), Tedges);
Tcexp = histc(Tol.Tmax(Imar&Itmc), Tedges);
Thab=Pobis.T_atPatm(:,(end-3):(end-1));
nThab = histc(Thab(:,2), Tedges);
for i=2:5
Taphi(:,i) = histc(Ta(:,i), Tedges1);
Ta2phi(:,i) = histc(Ta_combo(:,i), Tedges1);
end

Tol.refa=Tol.Phylum;
Iphy=[1 2 4 5 7:9]; % all animals
Iphy=[2 7 9]; % only chordate/mollust/crustac.
lref=unique(Tol.refa(Imar));

for i=1:length(lref)
    Iref=strcmp(lref{i},Tol.refa);
    frq_tc(:,i) = histc(Tol.Tmax(Iref),Tedges);
end
frq_tc=frq_tc(:,Iphy);
lref=lref(Iphy);

SST_day=ncread('~/Dropbox/DATA/sst.1998-2006.s04_day_pfv50_woagrid.nc','SST');
SSTc_day=calc_clim(SST_day);
SST_night=ncread('~/Dropbox/DATA/sst.1998-2006.s04_night_pfv50_woagrid.nc','SST');
SSTc_night=calc_clim(SST_night);
SST_mean=ncread('~/Dropbox/DATA/fv01_woagrid.nc','analysed_sst'); SST_mean=SST_mean-Tref;
SSTc_mean=calc_clim(SST_mean);

Tedges2=0:0.5:50; 
Tcent2=Tedges2(2:end)-mean(diff(Tedges2))/2;
ga=repmat(WOA.a3d(:,:,1),[1 1 12]);
for i=1:length(Tedges2)-1
    I=find(SSTc_day>Tedges2(i) & SSTc_day<Tedges2(i+1));
    fSST_day(i)=nansum(ga(I))./nansum(ga(:));
    I=find(SSTc_night>Tedges2(i) & SSTc_night<Tedges2(i+1));
    fSST_night(i)=nansum(ga(I))./nansum(ga(:));
    I=find(SSTc_mean>Tedges2(i) & SSTc_mean<Tedges2(i+1));
    fSST_mean(i)=nansum(ga(I))./nansum(ga(:));
end

I=find(WOA.z<=100);
Tsfc=WOA.temp(:,:,I,:);
vol=repmat(WOA.v3d(:,:,I),[1 1 1 12]);
for i=1:length(Tedges2)-1
    I=find(Tsfc>Tedges2(i) & Tsfc<Tedges2(i+1));
    fSST(i)=nansum(vol(I))./nansum(vol(:));
end

% plot

figure(5); clf reset;

py=[0.3 0.8];
subplot(321)
bar(Tedges1,Tcphi,'b'); hold on; axis([0 50 0 10])
plot(Tcent2,fSST/max(fSST)*max(Tcphi),'-','Color',[0.5 0.5 0.5])
plot(Tcent2,fSST_day/max(fSST_day)*max(Tcphi),':','Color',[0.5 0.5 0.5])
h3=legend('\Phi=1','Location','NorthWest');
xlabel('Temperature'); ylabel('# Species')

subplot(323) 
bar(Tedges,frq_tc,'stacked'); hold on; axis([0 50 0 250])
h1=legend(lref,'Location','NorthWest');set(h1,'FontSize',8)
plot(Tcent2,fSST/max(fSST)*max(sum(frq_tc,2)),'-','Color',[0.5 0.5 0.5])
plot(Tcent2,fSST_day/max(fSST_day)*max(sum(frq_tc,2)),':','Color',[0.5 0.5 0.5])
xlabel('Temperature'); ylabel('# Species')
[cmap]=cbrewer('seq','YlGn',5);
colormap(cmap)

subplot(325)
py=[0.3 0.8]*1.5;
h1=bar(Tedges1,Ta2phi(:,4),'FaceColor','r'); %'stacked');hold on;
hold on; axis([0 50 0 12])
h4=legend('\Phi=\Phi_{crit}','Location','NorthWest');
plot(Tcent2,fSST/max(fSST)*max(Ta2phi(:,4)),'-','Color',[0.5 0.5 0.5])
plot(Tcent2,fSST_day/max(fSST_day)*max(Ta2phi(:,4)),':','Color',[0.5 0.5 0.5])
xlabel('Temperature'); ylabel('# Species')

if doPrint
print -djpeg MI_traits_Fig5.jpg
print(gcf,'-depsc2', '-painters', 'MI_traits_Fig5.eps');
end
