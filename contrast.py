from periodictable import xsf
import periodictable as pt
import numpy as np
import math
import pretty_errors

def func(tinfrac,bkg = None):


    bd = {'dmf':['C3H7NO',.944], 'hexane': ['C6H8',.661]}
    

    indium = 1-tinfrac
    tin = tinfrac
    ox = 3/2*indium+2*tin

    total = np.sum([indium,tin,ox])

    in_frac = indium/total
    sn_frac = tin/total
    o_frac = ox/total

    unit = f'In{in_frac}Sn{sn_frac}O{o_frac}'
    
    
    print(unit)

    compound = pt.formula('In2O3')

    compound.density = 7.14

    nc_sld, mu = xsf.xray_sld(unit, density = 7.14, energy = 8.04)  #E10 / cm^2

    if bkg is not None:
        bkg_sld, mu = xsf.xray_sld(bd[bkg][0], density = bd[bkg][1], energy = 8.04) #E10/cm^2

        return((nc_sld - bkg_sld)**2)


if __name__ == '__main__':
    func(0.0, 'hexane')