function [ Phic, FratioT, FratioO] = get_Fstats( Pc4d )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


FvPhi=Pc4d.Fpz(Pc4d.kmax,:);
[Fmax,j]=max(FvPhi);
if j<length(Pc4d.prange) && j>1;
    Phic=Pc4d.prange(j);
else
    Phic=nan;
end
FratioT=max(Pc4d.FratioT);
%keyboard
FratioO=max(Pc4d.FratioO);

%max(nanmax(Pc4d.Fpz,[],2)./nanmax(Pc4d.Foz,[],2))

%zl=Pc4d.zl;
