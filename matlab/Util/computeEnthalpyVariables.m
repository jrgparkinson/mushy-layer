function [T,Sl, Ss, chi, HS, HE, HL] = computeEnthalpyVariables(H,S,p)

%Bounding energy


HS = p.cp*(p.thetaEutectic + max(0, -(p.C+S)/p.pc));
HE = (p.C + S)./(p.C + 1) * (p.St + p.thetaEutectic*(1-p.cp)) + p.cp*p.thetaEutectic;
HL = p.St - S + p.thetaEutectic + p.ThetaEutectic;

T = NaN*H; Sl = NaN*H; Ss = NaN*H; chi = NaN*H;

for i = 1:length(H)

    if H(i) <= HS(i)
       T(i) = H(i)/p.cp;
       Sl(i) = NaN;
       Ss(i) = S(i);
       chi(i) = 0;
       
    elseif H(i) <= HE(i)
         T(i) = p.thetaEutectic;
       Sl(i) = p.ThetaEutectic;
       
       chi(i) = (H(i)-p.thetaEutectic*p.cp)/(p.St + p.thetaEutectic*(1-p.cp));
       Ss(i) = S(i)/(1-chi(i));
    elseif H(i) < HL
        A = p.C*(p.cp-1) + p.St*(p.pc-1);
        B = H(i)*(1-p.pc) - p.St*p.pc + p.C*(1-2*p.cp) - S(i)*(p.cp-1);
        C = H(i)*p.pc + p.cp*(p.C + S(i));
        chi(i) = (-B-sqrt(B^2-4*A*C))/(2*A);
       
        Sl(i) = (S(i) + p.C*(1-chi(i)))/(chi(i) + (1-chi(i))*p.pc);
         T(i) = -Sl(i);
       
       Ss(i) = (p.pc*S(i) - p.C*(chi(i)))/(chi(i) + (1-chi(i))*p.pc);
       
    else
         T(i) = H(i) - p.St;
       Sl(i) = S(i);
       Ss(i) = NaN;
       chi(i) = 1;
    end
    
end


end