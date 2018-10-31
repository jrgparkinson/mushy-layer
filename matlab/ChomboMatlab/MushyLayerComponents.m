classdef MushyLayerComponents
    properties
        enthalpy
        enthalpySolid
        enthalpyEutectic
        enthalpyLiquid
        composition
        theta
        thetaLiquidus
        thetaSolidus
        porosity
        porosityEutectic
        compositionLiquid
        compositionSolid
        thetaForcing
        liquidCompositionGrad
        solidFraction
        steadyStateImbalance
        Tanalytic
        ThetaLAnalytic
        solidFractionTrue
        thetaTrue
        ThetaAnalytic
        enthalpyAnalytic
        enthalpyAdvection
        thetaLaplacian
        streamFunction
        resid
        concSource
        ThetaDiffusion
        ThetaDiffusionN
        thetaBcoef
        ThetaBCoef
        permeability
        pressure
        divU
        ThetaFrameAdvection
        ThetaSSource
        ThetaPorositySource
        pressureErr
        pressureAnalytic
        divUstarErr
        thetaForcingErr
        thetaForcingAnalytic
        CFInterpTest
        divUErr
        divUstar
        lambda
        lambdaPostCorr
        
        xFluidVel
        yFluidVel
        xFluidVelPred
        yFluidVelPred
        xGradP
        yGradP
        
        xFluidVelAnalytic
        yFluidVelAnalytic
        xFluidVelErr
        yFluidVelErr
        xGradPErr
        yGradPErr
        
    end
    methods
        function obj = MushyLayerComponents(comp_names)
%             obj.enthalpy = 1;
%             obj.enthalpySolid = 2;
%             obj.enthalpyEutectic = 3;
%             obj.enthalpyLiquid = 4;
%             obj.composition = 5;
%             obj.theta = 6;
%             obj.thetaLiquidus = 7;
%             obj.thetaSolidus = 8;
%             obj.porosity =9;
%             obj.porosityEutectic = 10;
%             obj.compositionLiquid = 11;
%             obj.compositionSolid = 12;
%             obj.thetaForcing = 13;
%             obj.liquidCompositionGrad = 14;
%             obj.solidFraction = 15;
%             obj.steadyStateImbalance = 16;
%             obj.Tanalytic = 17;
%             obj.ThetaLAnalytic = 18;
%             obj.solidFractionTrue = 19;
%             obj.thetaTrue = 20;
%             obj.ThetaAnalytic = 21;
%             obj.enthalpyAnalytic = 22;
%             obj.enthalpyAdvection = 23;
%             obj.thetaLaplacian = 24;
%             obj.streamFunction = 25;
%             obj.resid = 26;
%             obj.concSource = 27;
%             obj.ThetaDiffusion = 28;
%             obj.ThetaDiffusionN = 29;
%             obj.thetaBcoef = 30;
%             obj.ThetaBCoef = 31;
%             obj.permeability = 32;
%             obj.pressure = 33;
%             obj.divU = 34;
%             obj.ThetaFrameAdvection = 35;
%             obj.ThetaSSource = 36;
%             obj.ThetaPorositySource = 37;
%             obj.pressureErr = 38;
%             obj.pressureAnalytic = 39;
%             obj.divUErr = 40;
%             obj.thetaForcingErr = 41;
%             obj.thetaForcingAnalytic=42;
%             obj.CFInterpTest = 43;
%             obj.divUstar = 44;
%             obj.divUstarErr = 45;
%             obj.lambda = 46;
%             obj.lambdaPostCorr = 47;
%             
%             
%             vectorVarsInit = 48; % This should be one greater than the largest scalar var index
%             
%             obj.xFluidVel = vectorVarsInit;
%             obj.zFluidVel = vectorVarsInit+1;
%             obj.xFluidVelPred = vectorVarsInit+2;
%             obj.zFluidVelPred = vectorVarsInit+3;
%             obj.xGradP = vectorVarsInit+4;
%             obj.zGradP = vectorVarsInit+5;
%             
%             obj.xFluidVelAnalytic = vectorVarsInit+6;
%             obj.zFluidVelAnalytic = vectorVarsInit+7;
%             obj.xFluidVelErr = vectorVarsInit+8;
%             obj.zFluidVelErr = vectorVarsInit+9;
%             obj.xGradPErr = vectorVarsInit+10;
%             obj.zGradPErr = vectorVarsInit+11;
            
            if nargin > 0
                obj = obj.parseCompList(comp_names);
            end
            
        end
        
        function obj = parseCompList(obj, comp_names)
        
            obj.enthalpy = find(strcmp(comp_names,'Enthalpy')); 
            obj.enthalpySolid = find(strcmp(comp_names,'Enthalpy solidus'));
            obj.enthalpyEutectic = find(strcmp(comp_names,'Enthalpy eutectic'));
            obj.enthalpyLiquid = find(strcmp(comp_names,'Enthalpy liquidus'));
            obj.composition = find(strcmp(comp_names,'Composition'));
            obj.theta = find(strcmp(comp_names,'theta'));
            obj.thetaLiquidus = find(strcmp(comp_names,'theta liquidus'));
            obj.thetaSolidus = find(strcmp(comp_names,'theta solidus'));
            obj.porosity = find(strcmp(comp_names,'Porosity'));
            obj.porosityEutectic = find(strcmp(comp_names,'Eutectic Porosity'));
            obj.compositionLiquid = find(strcmp(comp_names,'Liquid Composition'));
            obj.compositionSolid = find(strcmp(comp_names,'Solid Composition'));
            obj.thetaForcing = find(strcmp(comp_names,'theta forcing'));
            obj.liquidCompositionGrad = find(strcmp(comp_names,'Liquid Composition Gradient'));
            obj.solidFraction = find(strcmp(comp_names,'Solid Fraction'));
            obj.steadyStateImbalance = find(strcmp(comp_names,'Steady state imbalance'));
            obj.Tanalytic = find(strcmp(comp_names,'Analytic Temperature'));
            obj.ThetaLAnalytic = find(strcmp(comp_names,'Analytic liquid concentration'));
            obj.solidFractionTrue = find(strcmp(comp_names,'Analytic Solid Fraction'));
            obj.thetaTrue = find(strcmp(comp_names,'Analytic theta'));
            obj.ThetaAnalytic = find(strcmp(comp_names,'Theta analytic'));
            obj.enthalpyAnalytic =find(strcmp(comp_names,'Enthalpy analytic'));
            obj.enthalpyAdvection = find(strcmp(comp_names,'Enthalpy advection'));
            obj.thetaLaplacian = find(strcmp(comp_names,'theta laplacian'));
            obj.streamFunction = find(strcmp(comp_names,'Stream function'));
            obj.resid = find(strcmp(comp_names,'Residual'));
            obj.concSource = find(strcmp(comp_names,'Concentration Eqn Source'));
            obj.ThetaDiffusion = find(strcmp(comp_names,'Theta Diffusion at n'));
            obj.ThetaDiffusionN = find(strcmp(comp_names,'Theta Diffusion'));
            obj.thetaBcoef = find(strcmp(comp_names,'theta B coeff'));
            obj.ThetaBCoef = find(strcmp(comp_names,'Theta B coeff'));
            obj.permeability = find(strcmp(comp_names,'Permeability'));
            obj.pressure = find(strcmp(comp_names,'Pressure'));
            obj.divU = find(strcmp(comp_names,'Divergence U'));
            obj.ThetaFrameAdvection = find(strcmp(comp_names,'Theta frame adv'));
            obj.ThetaSSource = find(strcmp(comp_names,'Theta solid source'));
            obj.ThetaPorositySource = find(strcmp(comp_names,'Theta porosity source'));
            obj.pressureErr = find(strcmp(comp_names,'Pressure Error'));
            obj.pressureAnalytic = find(strcmp(comp_names,'Pressure analytic'));
            obj.divUErr = find(strcmp(comp_names,'Div U rror'));
            obj.thetaForcingErr = find(strcmp(comp_names,'theta forcing error'));
            obj.thetaForcingAnalytic= find(strcmp(comp_names,'theta forcing analytic'));
            obj.CFInterpTest = find(strcmp(comp_names,'CF interp test'));
            obj.divUstar = find(strcmp(comp_names,'div U star'));
            obj.divUstarErr = find(strcmp(comp_names,'div U star err'));
            obj.lambda = find(strcmp(comp_names,'lambda passive tracer'));
            obj.lambdaPostCorr = find(strcmp(comp_names,'lambda after correction'));
            
            obj.xFluidVel = find(strcmp(comp_names,'xFluid velocity'));
            obj.yFluidVel = find(strcmp(comp_names,'yFluid velocity'));
            obj.xFluidVelPred = find(strcmp(comp_names,'xFluid velocity prediction'));
            obj.yFluidVelPred = find(strcmp(comp_names,'yFluid velocity prediction'));
            obj.xGradP = find(strcmp(comp_names,'xGrad pressure'));
            obj.yGradP = find(strcmp(comp_names,'yGrad pressure'));
            obj.xFluidVelAnalytic = find(strcmp(comp_names,'xFluid vel analytic'));
            obj.yFluidVelAnalytic = find(strcmp(comp_names,'yFluid vel analytic'));
            obj.xFluidVelErr = find(strcmp(comp_names,'xFluid vel error'));
            obj.yFluidVelErr = find(strcmp(comp_names,'yFluid vel error'));
            obj.xGradPErr = find(strcmp(comp_names,'xGrad Pressure Err'));
            obj.yGradPErr = find(strcmp(comp_names,'zGrad Pressure Err'));
        end
    end
end