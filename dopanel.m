
spec1=Par.Species_unq(i);
fname=strrep(spec1{1},' ','_');
[Pdat] = extract_Spec_MIparm(Par,spec1);

load(['output_OBIS/output_' fname '_mod'])
v2struct(Pc4d)
v2struct(Pcssp)

TObin_mod=log10(TObin_gt1);
TObin_mod(isnan(bin.P))=-1;

% Simple Met.Idx.
pcrit1 = @(b,x)b(1).*[exp(-b(2)*x(:,1))];

pcolor(Tcent,Ocent,TObin_mod');shading flat;hold on;colorbar;
xlabel('Temperature');ylabel('pO2 [atm]')
plot(Par.P.temp(i,:),Par.P.pcrit(i,:)/atm2kPa,'k*');
TT=(1/kb)*(1./(Tcent+T0)-1./(Tref+T0));
pfitEx1=pcrit1([1/Par.P.est(i,1) Par.P.est(i,3)],[TT']); 
plot(Tcent,pfitEx1,'k--');
plot(Tcent,pfitEx1*pc_F3,'k--');
hold on;
xlim(xl)
caxis(cax)
title([spec1]);
