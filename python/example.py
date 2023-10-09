import icegiant as ig
import numpy as np
import astropy.constants as cons
import astropy.units as unit

#instantiate the parameter object
parameters = ig.params()


#populate the parameter container 
parameters.mission ="LISA"
parameters.thetaL = np.pi/2
parameters.phiL = 0.4
parameters.thetaS = np.pi/3
parameters.phiS = 0.5 
parameters.M1 = (0.2*unit.M_sun).to(unit.kg).value #mass of first WD
parameters.M2 = parameters.M1 #mass of second WD
parameters.MP = (1.*unit.M_jup).to(unit.kg).value #mass of circumbinary planet
parameters.Tobs = (4*unit.yr).to(unit.s).value #observation time
parameters.Larm = 2.5e9 #length of LISA arm
parameters.NC = 2 ##number of LISA channels
parameters.thetaP = np.pi/2
parameters.phiP = np.pi/2
parameters.freqGW = 10e-3
parameters.sourceDistance = (1e3*unit.pc).to(unit.m).value #distance to source
parameters.ig_direction = 0
parameters.LISAAlpha = 2
parameters.lightTwoWayTime = 15000 #time for light to travel to spacecraft and back
parameters.DerivativeDelta = 1e-6 #delta X for calculating derivatives



#instantiate Binary 
bi = ig.Binary(parameters)

print(bi.strain(parameters,3.))

