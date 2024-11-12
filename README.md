# MassfitCorrectionsZ
Muon momentum bias corrections from Z mass fits

run_massloop_data.py runs all the steps in the correct order:
./massscales_data
./massfit
./resolfit 
  -> once 
./massscales_data (command updated to use results of the previous iteration)
./massfit 
./resolfit 
  -> several iterations

massscales_data.cpp fits the Z mass in 4D bins (muon eta+, pt+, eta-, pt-) for data and MC extracting the mass scale and resolution biases. It can be used iteratively with massfit.cpp and resolfit.cpp to correct muon momenta in MC to get mass distributions closer to the data.

massfit.cpp can be ran in:
  data mode -> takes the mass scale biases per 4D bin and fits for the pT scale correction parameters A,e,M per eta bin 
  OR toys mode -> generates mass scale biases from dummy AeM biases and fits for AeM from them (a closure test)

resolfit.cpp can be ran in:
  data mode -> takes the mass width biases per 4D bin and fits for the pT resolution correction parameters c,d per eta bin
  OR toys mode -> generates mass width biases from dummy cd biases and fits for cd from them (a closure test)

