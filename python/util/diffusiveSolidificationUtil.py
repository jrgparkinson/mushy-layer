# Assumption through this file is that 
# at z=1, theta = theta_eutectic
# with the domain extending down to z=-infinity, where theta = theta_infinity
# params should be a dict containing the fields:
#   nonDimVel, compositionRatio, stefan
# Calculations follow Worster (1991)

import math

# Compute the mushy layer height
def computeMushyH(params, thetaInf, thetaInterface, zEutectic=1.0):
    CR = float(params['parameters.compositionRatio'])
    St = float(params['parameters.stefan'])
    V  = float(params['parameters.nonDimVel'])
    
    A = 0.5*(CR + St + thetaInf)
    B = math.sqrt(A**2 - CR *thetaInf - St)
    
    alpha = A+B
    beta = A-B
    
    hV = V*zEutectic - ((alpha-CR)/(alpha-beta)) * math.log(alpha/(alpha-thetaInterface)) - ((CR-beta)/(alpha-beta))*math.log(beta/(beta-thetaInterface))
    h = hV/V
    
    return h

# Compute far field temperature
def compute_thetaInf(params, HBottom, zbottom, thetaInterface):
    CR = float(params['parameters.compositionRatio'])
    St = float(params['parameters.stefan'])
    V  = float(params['parameters.nonDimVel'])
    
    thetaBottom = HBottom - St
    
    thetaIncrement = 1e-4
    
    thetaEutectic = thetaInterface - 1.0
    
    thetaInf = thetaBottom # First guess
   
    
    while thetaInf < thetaBottom + 10.0:
    
        zEutectic_predicted = compute_eutecticPosition(params, thetaInf, thetaInterface, thetaEutectic)
        h_predicted = computeMushyH(params, thetaInf, thetaInterface, zEutectic_predicted)
        thetaBottom_predicted = compute_thetaBottom(params, thetaInf, thetaInterface, zbottom, V, h_predicted)
        
        print('zEutectic predicted: ' + str(zEutectic_predicted) + ', h predicted: ' + str(h_predicted) + ', thetaBottom predicted: ' + str(thetaBottom_predicted))
        
        if abs(thetaBottom_predicted - thetaBottom) < thetaIncrement*10:
            break

        thetaInf = thetaInf + thetaIncrement
        
    print('Determined theta infinity = ' + str(thetaInf))
        
    return thetaInf

# Compute temperature at the bottom of the domain
def compute_thetaBottom(params, thetaInf, thetaInterface, zbottom, V, h=0):
    if h == 0:
        thetaEutectic = thetaInterface - 1.0
        zEutectic = compute_eutecticPosition(params, thetaInf, thetaInterface, thetaEutectic)
        h = computeMushyH(params, thetaInf, thetaInterface, zEutectic)
        
    return thetaInf + (thetaInterface - thetaInf)*math.exp(V*(zbottom-h))
    

# Compute the enthalpy at the bottom of the domain given the far field temperature
def compute_HBottom(params, thetaInf, zbottom, thetaInterface):
    St = float(params['parameters.stefan'])
    V  = float(params['parameters.nonDimVel'])
    
    theta_bottom = compute_thetaBottom(params, thetaInf, thetaInterface, zbottom, V)
    
    HBottom = theta_bottom + St
    
    return HBottom
    
    
def compute_eutecticPosition(params, thetaInf, thetaInterface, thetaEutectic):
    V  = float(params['parameters.nonDimVel'])
    St = float(params['parameters.stefan'])
    CR = float(params['parameters.compositionRatio'])
    thetaTop = float(params['parameters.topEnthalpy'])
    
    A = 0.5*(CR + St + thetaInf)
    dThetadz_eutectic = (thetaEutectic**2 - 2*thetaEutectic*A+CR*thetaInf + St)/(V*(CR-thetaEutectic))
    
    #print('dThetadz_eutectic: ' + str(dThetadz_eutectic))
    
    const = V*thetaEutectic + dThetadz_eutectic
    
    #print('Const: ' + str(const))
        
    zEutectic = 1 - math.log((V*thetaTop - const)/(V*thetaEutectic - const))/V

    return zEutectic
