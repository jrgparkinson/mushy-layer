function [ averagedTdH, averagedSldS ] = calcPhaseDiagDerivs( St, CR, cp )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
theta_eutectic = -1;
    Theta_eutectic = 1;
    
    pc = 1e-5;
    
    H_lower = 1.0;
    H_warm = 5.0;
    
    S_lower = 0;
    S_max = Theta_eutectic;
    
    
    
    Hvec = linspace(H_lower, H_warm, 50);
    Svec = linspace(S_lower, S_max, 50);
    
    [S, H] = meshgrid(Svec, Hvec);
    
    T = NaN*H;
    
    Sl = T; chi = T; Ss = T;
    TLinear = T; chiLinear = T; SlLinear = T; SsLinear = T;
    dTdHMush = T; dSldSMush = T; dChidHMush = T; dChidSMush = T;
    
    % Increased salinity = increased porosity
    chi_eutectic = (CR+S)/(CR+Theta_eutectic);
    
    % Bounding energy
    H_s = cp*(theta_eutectic + max(0, -(S+CR)/pc)) ; %H_s = -1
    H_e = chi_eutectic*(St+theta_eutectic*(1-cp)) + cp*theta_eutectic;
    H_l = St - S + theta_eutectic + Theta_eutectic;
    
    % Compute T, Sl etc.
    for H_i = 1:length(Hvec)
        for S_i = 1:length(Svec)
            
            if Hvec(H_i) <= H_s(H_i, S_i)
                T(H_i, S_i) = Hvec(H_i)/cp;
                Sl(H_i, S_i) = 0;
                chi(H_i, S_i) = 0;
                
            elseif Hvec(H_i) <= H_e(H_i, S_i)
                T(H_i, S_i) = theta_eutectic;
                Sl(H_i, S_i) = Theta_eutectic;
                
                chi(H_i, S_i) = (Hvec(H_i) - theta_eutectic*cp)/(St + theta_eutectic*(1-cp));
                
            elseif  Hvec(H_i) < H_l(H_i, S_i)
                A = CR*(cp-1) + St*(pc-1);
                B = Hvec(H_i) *(1-pc) - St*pc + CR*(1-2*cp) - Svec(S_i)*(cp-1);
                C = Hvec(H_i)*pc + cp*(CR + Svec(S_i));
                
                chi(H_i, S_i) = (-B - sqrt(B^2 - 4*A*C))/(2*A);
                Sl(H_i, S_i) = (S(H_i, S_i) + CR*(1-chi(H_i, S_i)))/ ...
                    (chi(H_i, S_i) + (1-chi(H_i, S_i))*pc);
                
                T(H_i, S_i) = - Sl(H_i, S_i);
                
                
                dChidHMush(H_i, S_i) = -(1/(2*A))*(1+B/sqrt(B^2 - 4*A*C));
                dChidSMush(H_i, S_i) = -(1/(2*A))*((1-cp) + (B*(1-cp)-2*A*cp)/(sqrt(B^2-4*A*C)));
                
                dTdHMush(H_i, S_i) = (1/chi(H_i, S_i)^2)*dChidHMush(H_i, S_i)* ...
                    (Svec(S_i) + CR);
                dSldSMush(H_i, S_i) = -(1/chi(H_i, S_i)^2)*dChidSMush(H_i, S_i)* ...
                    (Svec(S_i) + CR) + 1/chi(H_i, S_i);
                
                
            else
                T(H_i, S_i) =  Hvec(H_i) - St;
                Sl(H_i, S_i) = Svec(S_i);
                chi(H_i, S_i) = 1;
            end
        end
    end
    
    averagedTdH = mean(mean(dTdHMush(~isnan(dTdHMush))));
    averagedSldS = mean(mean(dSldSMush(~isnan(dSldSMush))));
    

end

