import hapi as hp
import matplotlib.pyplot as plt

hp.db_begin("data")
#hp.fetch("H2O",1,1,3400,4100)
hp.select('H2O',ParameterNames=('nu','sw'),Conditions=('between','nu',4000,4100))


nu,coef =  hp.absorptionCoefficient_Voigt(((1,1),),'H2O', OmegaStep=1,HITRAN_units=False,GammaL='gamma_self', Environment={'p':1,'T':296.})

plt.plot(nu,coef)
plt.show()