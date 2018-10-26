function [heatAdvection, heatDiffusion, latentHeat, TFrameAdvection, ...
   saltAdvection, saltDiffusion, liquidSalinityFrame, solidSalinityFrame, ...
  vorticityDiffusion, baroclinicTorque, vorticityPermeability] = ...
  computeFields(finalPlotFile, frameAdv, St, CR, Ra, Le)

% Always apply this check
    if length(finalPlotFile.levelArray()) == 0
        return
    end

   

    T = finalPlotFile.dataForComp(finalPlotFile.components.Temperature).';
    chi = finalPlotFile.dataForComp(finalPlotFile.components.Porosity).';
    H = finalPlotFile.dataForComp(finalPlotFile.components.Enthalpy).';
    S = finalPlotFile.dataForComp(finalPlotFile.components.Bulkconcentration).';
    Sl = finalPlotFile.dataForComp(finalPlotFile.components.Liquidconcentration).';
    
    U = finalPlotFile.dataForComp(finalPlotFile.components.xAdvectionvelocity).';
    V = finalPlotFile.dataForComp(finalPlotFile.components.yAdvectionvelocity).';
    
      
      
    psi = finalPlotFile.getStreamfunction(1e4);
      
    %psi2 =  finalPlotFile.getStreamfunction(1e5);
    %psiDiff = (psi2-psi)./psi;
    % psiMaxFracDiff = max(max(psiDiff));
     chiSl = chi.*Sl;
     
    perm = chi.^3;
    
    
  
    
     % Do some conversion to Wells Units
   
    Sl = Sl + 1;
    S = S + 1;
    
    H = H - 1;
    T = T - 1;
    
     [X, Y] = finalPlotFile.grid();
    dx = X(1, 2) - X(1, 1);
    x = X(1, :);
    y = Y(:, 1);
    
    % Heat equation
    [dTdx, dTdz] = gradient(T, dx);
    [dHdx, dHdz] = gradient(H, dx);
    [dchidx, dchidz] = gradient(chi, dx);
    
    heatAdvection = U.*dTdx + V.*dTdz;
    heatDiffusion = 4*del2(T, dx);
    latentHeat = frameAdv*St*dchidz;
    TFrameAdvection = frameAdv*dTdz;
    
    
    % Salt equation
    [dSldx, dSldz] = gradient(Sl, dx);
    [dSdx, dSdz] = gradient(S, dx);
    
    saltAdvection = U.*dSldx + V.*dSldz;
    saltDiffusion = divergence(chi.*dSldx, chi.*dSldz)/Le;
    [~, liquidSalinityFrame] = gradient(chi.*Sl * frameAdv, dx);
    solidSalinityFrame = - dchidz * frameAdv * CR;
    
    % Momentum equation
    [dpsidx, dpsidz] = gradient(psi, dx);
    [dpermdx, dpermdz] = gradient(perm, dx);
    vorticityDiffusion = 4*del2(psi, dx);
    baroclinicTorque = Ra*dTdx.*perm;
    vorticityPermeability = (dpsidx.*dpermdx + dpsidz.*dpermdz)./perm;
    
   
%     figure();
%     m = 2; n = 2;
%     
%     subplot(m,n, 1);
%     h = pcolor(perm);  set(h, 'EdgeColor', 'none');
%     colorbar(); ylabel('Pi');
%     
%     subplot(m,n, 2);
%     h = pcolor(dpermdx);  set(h, 'EdgeColor', 'none');
%     colorbar(); ylabel('d(Pi)/dx');
%    
%     subplot(m,n, 3);
%     %h = pcolor(T);  set(h, 'EdgeColor', 'none');
%     contour(T);
%     colorbar(); ylabel('T');
%     
%     subplot(m,n, 4);
%     h = pcolor(dTdx);  set(h, 'EdgeColor', 'none');
%     colorbar(); ylabel('dT/dx');
    
%     Uerr = U + dpsidz;
%     Verr = V + dpsidx;
%     
%     figure();
%     m = 2;
%     n = 1;
    
%     subplot(m, n, 1);
%     h = pcolor(Uerr);  set(h, 'EdgeColor', 'none');
%     ylabel('U err');
%     colorbar();
%     
%     subplot(m, n, 2);
%     h =pcolor(Verr);  set(h, 'EdgeColor', 'none');
%     ylabel('V err');
%     colorbar();
    
    

    % Let's try and calculate these things correctly
%     for i=2:length(x)-1
%         for j=2:length(y)-1
%             
%             % Remember that x is the second index
%             
%             gradT_x = (T(j, i+1) - T(j, i-1))/(2*dx);
%             gradT_y = (T(j+1, i) - T(j-1, i))/(2*dx);
%             gradChi_y = (chi(j+1, i) - chi(j-1, i))/(2*dx);
%             
%             gradSl_x = (Sl(j, i+1) - Sl(j, i-1))/(2*dx);
%             gradSl_y = (Sl(j+1, i) - Sl(j-1, i))/(2*dx);
%             
%         
%             u = U(j, i); v = V(j,i);
%             
%             heatAdvection(j, i) = u*gradT_x + v*gradT_y;
%             heatDiffusion(j, i) = (T(j, i+1) + T(j, i-1) + T(j-1, i) + T(j+1, i) - 4*T(j, i))/(dx*dx);
%             latentHeat(j, i) = frameAdv*St*gradChi_y;
%             TFrameAdvection(j, i) = frameAdv*gradT_y;
%             
%             
%             
%             
%             saltAdvection(j, i) = -heatAdvection(j, i);
%             
%             %saltDiffusion(j, i) = 0; % This is really small anyway
%             
%             liquidSalinityFrame(j, i) = frameAdv*(chiSl(j+1,i) - chiSl(j-1,i))/(2*dx);
%             solidSalinityFrame(j, i) = frameAdv * CR * frameAdv*(chi(j+1,i) - chi(j-1,i))/(2*dx);
%             
%             vorticityDiffusion(j, i) = (psi(j, i+1) + psi(j, i-1) + psi(j-1, i) + psi(j+1, i) - 4*psi(j, i))/(dx*dx);
%             baroclinicTorque(j, i) = Ra*dTdx.*perm;
%             vorticityPermeability(j, i) = (dpsidx.*dpermdx + dpsidz.*dpermdz)./perm;
%             
%         end
%     end
    



end