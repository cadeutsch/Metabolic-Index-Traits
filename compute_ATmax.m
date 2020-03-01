function [ Tc, Ta, Ta2 ] = compute_ATmax(Par,Pobis)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global kb Rgas T0 Tref Patm

Ac=Par.P.est(:,1); 
Ec=Par.P.est(:,3); 
dEdT=Par.dEdT; %0.6/30;

% x=Par.tmedP; y=Par.P.est(:,3);
% [B,BINT,R,RINT,STATS] = regress(y,[x x./x]);
% dEdT=B(1);

myfun = @(x,c) c(1)*exp(c(2) * (c(3) + c(4)*(x-c(5))) * (1/x - 1/c(6))) -1;

clear AT*
for i=1:length(Ec)
    if Ec(i)>0.01
        
        c0 = [Ac(i)*Patm 1/kb Ec(i) 0 Par.tmedP(i)+T0 Tref+T0];
        fun0 = @(x) myfun(x,c0);    % function of x alone
        sol0 = fzero(fun0,T0+30);
        
        c1 = [Ac(i)*Patm 1/kb Ec(i) dEdT(i) Par.tmedP(i)+T0 Tref+T0];
        fun1 = @(x) myfun(x,c1);    % function of x alone
        sol1 = fzero(fun1,T0+30);
        
        ATmax_phi1_Econst(i)=sol0-T0;
        ATmax_phi1_dEdT(i)=sol1-T0;
        
        if isfinite(Pobis.phi_crit(i))
            
            c0 = [Ac(i)*Patm/Pobis.phi_crit(i) 1/kb Ec(i) dEdT(i) Par.tmedP(i)+T0 Tref+T0];
            fun0 = @(x) myfun(x,c0);    % function of x alone
            sol0 = fzero(fun0,T0+30);
            
            c1 = [Ac(i)*Patm/nanmean(Pobis.phi_crit) 1/kb Ec(i) dEdT(i) Par.tmedP(i)+T0 Tref+T0];
            fun1 = @(x) myfun(x,c1);    % function of x alone
            sol1 = fzero(fun1,T0+30);
            
            ATmax_phic_Sp(i)=sol0-T0;
            ATmax_phic_mnSp(i)=sol1-T0;
        end
        
    end    
    
%     if i==293
%         keyboard
%     end
    
end

Tc=ATmax_phi1_dEdT; Tc(Tc==0)=nan; % resting ATmax
Ta=ATmax_phic_Sp; Ta(Ta==0)=nan; % active ATmax w/ species-specific Phi_crit
Ta2=ATmax_phic_mnSp; Ta2(Ta2==0)=nan; % active ATmax w/ species-mean Phi_crit


end

