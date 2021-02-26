prm = {
    'gammaE': 0.1, # E lorentzian width
    'gammaI': 0.1, # I lorentzian width
    'omegaE_bar': 1.5, # E mean angular velocity
    'omegaI_bar': 0.5, # I mean angular velocity

    'KEE': 0, # coupling strength
    'KEI': 0.5, 
    'KIE': 0.5,
    'KII': 0,
    'N': 2000, # amount of oscillators per species
    'r': 1,

    'D': 0, # noise strength

    'dt': 0.0001, # time step
    'tmax': 100, # max. simulation time
    'step': 67, # data recording interval

    'use_cpp': False
}

storage = 'temp'
