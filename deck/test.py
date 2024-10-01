import opm_CoSTA
import numpy as np

blackoil = opm_CoSTA.BlackOil('TRACER')
print(blackoil.ndof)

nocorr = np.zeros(blackoil.ndof)
nsol = np.zeros(blackoil.ndof)

for i in range(0,30):
  sol = blackoil.predict(nsol)
  nsol = blackoil.correct({}, sol, nocorr)

#sol = blackoil.predict(nsol)
#residual = blackoil.residual({}, sol, uprev)
#uprev[2:-1:3] = 1e-6
#nsol = blackoil.correct({}, sol, uprev)
