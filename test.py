import opm_CoSTA

#blackoil = opm_CoSTA.BlackOil('/home/akva/kode/opm/opm-tests/spe1/SPE1CASE1')
#blackoil.predict([])

one_phase = opm_CoSTA.OnePhase('/home/akva/kode/opm/opm-tests/spe1/SPE1CASE1_WATER')
one_phase.predict([])
