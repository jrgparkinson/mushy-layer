
def computePorosityMush(CR, T, S):
    Sl = -T
    chi = (S+CR)/(Sl+CR)
    return chi

def get_sea_ice_material_properties():

    properties = {'Se':230,  # eutectic salinity, g/kg
    'Te': -23 , # eutectic temperature, celcius
'kappa_l' : 1.25e-7,  # heat diffusivity m^2 s^-1
'eta' :1.55e-3 , # liquid viscosity kg m^-1 s^-1
'rho_l' : 1028  ,# liquid density kg m^-3
'beta' : 7.86e-4 , # haline contraction (g/kg)^-1
'g' : 9.8  ,# gravitational acceleration m/s
'alpha' : 3.87e-5,  # thermal expansion celcius^-1
                  'L': 333.4e3, # latent heat of freezing J / kg
                  'kl': 0.523, # heat conductivity (liquid) W/m/celcius
                  'ks': 2.22, # heat conductivity (solid) W/m/celcius
                  'cpl': 4185, # heat capacity (liquid) J/kg/celcius
                  'cps': 2212, # heat capacity (solid) J/kg/celcius
                  'liquidusSlope': -0.1 # llinearised iquidus slope (celcius / g/kg)
    }


    return properties


# Ttop = -5.0 # top temperature, celcius
# Tbottom = -2.9 # bottom temperature, celcius
# Si = 30.0 # initial liquid salinity, g/kg
# h = 1.0 # length scale, metres
# d = 1e-4 # hele-shaw cell gap size (0.1 mm)
# K0 = 1e-10 # reference permeability. This is maybe a bit small. c.f. Rees Jones and Worster GRL paper
def set_params(params, Ttop, Tbottom, S_top = 0.0, Si=30.0, h=1.0, d=1e-4, K0=1e-10,
#                Se=230,  # eutectic salinity, g/kg
#     Te = -23 , # eutectic temperature, celcius
# kappa_l = 1.25e-7,  # heat diffusivity m^2 s^-1
# eta = 1.55e-3 , # liquid viscosity kg m^-1 s^-1
# rho_l = 1028  ,# liquid density kg m^-3
# beta = 7.86e-4 , # haline contraction (g/kg)^-1
# g = 9.8  ,# gravitational acceleration m/s
# alpha = 3.87e-5  , # thermal expansion celcius^-1
               darcy_brinkman=False,
               properties=get_sea_ice_material_properties(),
               dim=2, periodic=True):
    # Physical constants

    params['dimensional_values'] = '"Ttop=%.3g celcius, Tbottom=%.3g celcius, Si=%.3g g/kg, ' \
                                   'L=%.3g metres, d=%.3g metres, K0=%.3g m^2"' % (Ttop, Tbottom, Si, h, d, K0)

    ##############################
    # Compute quantities derived from our dimensional inputs
    ###############################

    Cref = properties['Se']
    Ti = properties['liquidusSlope'] * Si # liquidus temperature for initial salinity
    delta_c = properties['Se'] - Si
    delta_T = Ti - properties['Te']

    ###############################
    # Compute equivalent dimensionless values
    ##############################

    thetaTop = (Ttop - properties['Te'])/delta_T
    thetaBottom = (Tbottom-properties['Te'])/delta_T
    Thetai = (Si - Cref)/delta_c
    CR = Cref/delta_c

    kappa_l = properties['kl'] / (properties['rho_l']*properties['cpl'])

    RaS = properties['beta']*properties['rho_l']*properties['g']*delta_c*h**3/(kappa_l * properties['eta'])

    Da = K0/h**2
    RmS = Da*RaS
    Le = 200
    Pr = properties['eta']/(properties['rho_l']*kappa_l)
    # timescale = h**2/kappa_l
    St = properties['L'] / (properties['cpl'] * abs(delta_T))
    ThetaTop = (S_top - properties['Se'])/delta_c
    # ThetaTop = -CR
    # ThetaTop = Thetai
    ThetaBottom = Thetai
    chiTop = computePorosityMush(CR, thetaTop, ThetaTop)
    HTop = St*chiTop + thetaTop
    HBottom = St + thetaBottom

    # ThetaTop = thetaTop*chiTop + CR* (1-chiTop)

    RaT = properties['alpha']*properties['rho_l']*properties['g']*delta_T*h**3/(kappa_l * properties['eta'])
    RmT = RaT*Da

    heleShawPerm = d**2/(12*K0)

    if darcy_brinkman:
        params['parameters.rayleighComp'] = RaS
        params['parameters.rayleighTemp'] = RaT
        params['parameters.darcy'] = Da
        params['parameters.prandtl'] = Pr

    else:

        params['parameters.rayleighComp'] = RmS
        params['parameters.rayleighTemp'] = RmT
        params['parameters.darcy'] = 0
        params['parameters.prandtl'] = 0


        # params['parameters.nonDimReluctance'] = reluctance
    params['parameters.heleShawPermeability'] = heleShawPerm
    params['parameters.stefan'] = St
    params['parameters.compositionRatio'] = CR
    params['parameters.lewis'] = Le

    params['parameters.heatConductivityRatio'] = properties['ks'] / properties['kl']
    params['parameters.specificHeatRatio'] = properties['cps']  / properties['cpl']

    # No flux at the top boundary breaks the projection (we can't get a divergence free field)
    #params['bc.bulkConcentrationHi'] = '1 1'  # no flux at top boundary
    #params['bc.bulkConcentrationHiVal']= '0 0'
    #params['bc.bulkConcentrationLoVal']= '0 -1'  # fixed value at bottom boundary

    bc_accuracy = 4 # num significant figures
    if dim == 2:
        bc_str = '0 %.' + str(bc_accuracy) +'g'  # 0 in x dir, some val in vertical dir
    else:
        bc_str = '0 0 %.' + str(bc_accuracy) +'g'

    if periodic:
        params['main.periodic_bc'] = '1 ' * (dim-1) + '0'
    else:
        params['main.periodic_bc'] = '0 ' * dim

    params['bc.bulkConcentrationHi'] = '1 ' * dim # no salt flux at either boundary
    params['bc.bulkConcentrationHiVal']= '0 ' * dim # no salt flux at either boundary
    params['bc.bulkConcentrationLoVal']= bc_str % ThetaBottom  # fixed value (ocean salinity) at bottom boundary

    params['bc.enthalpyHiVal']= bc_str % HTop
    params['bc.enthalpyLoVal']= bc_str % HBottom

    params['bc.temperatureHiVal']= bc_str % thetaTop
    params['bc.temperatureLoVal']= bc_str % thetaBottom

    params['bc.scalarHi']= '1 ' * (dim-1) + '0'
    params['bc.scalarLo']= '1 ' * (dim-1) + '2'

    params['bc.velHi'] = '0 ' * dim
    params['bc.velLo'] = '0 ' * dim

    if darcy_brinkman and params['parameters.heleShawPermeability'] * params['parameters.darcy'] < 10.0:
        print('Warning - for Darcy Brinkman need \\Pi_H * Da \ll 1')

    return params

if __name__ == "__main__":

    #
    dim = 3
    periodic = True

    material_properties = get_sea_ice_material_properties()
    Si = 30.0  # initial salinity (g/kg)

    Ttop = -10 # top temperature (celcius)
    initial_freezing_point = material_properties['liquidusSlope'] * Si
    Tbottom = initial_freezing_point + 0.3  # ocean temperature - initial freezing point  (celcius)

    h = 1.0  # box depth (m)
    d = 1e-4  # Hele-Shaw gap width (m)
    K0 = 1e-10  # Sea ice permeability

    darcy_brinkman = False

    params = {}
    p = set_params(params, Ttop, Tbottom, Si=Si, h=h, d=d, K0=K0,
                   darcy_brinkman=darcy_brinkman, properties = material_properties, dim=dim,
                   periodic=periodic)

    for key in p:
        try:
            val = '%.3g' % p[key]
        except TypeError:
            val = '%s' % p[key]


        print('%s=%s' % (key, val))
