function [ Spec ] = extract_Spec_MIparm( Preg, sp_name )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
      
ss=strcmp(Preg.Species_unq,sp_name);
    Spec.sp_name = sp_name;

try
    % MI parameters
    Spec.Ao = Preg.P.est(ss,1);
    Spec.eps = Preg.P.est(ss,2);
    Spec.Eo = Preg.P.est(ss,3);
    Spec.B = Preg.B(ss);
    Spec.Beco = Preg.Beco(ss);
    if isnan(Spec.B) || Spec.B==0
        Spec.B=1;
    end
catch
    disp('Species not found')
end
