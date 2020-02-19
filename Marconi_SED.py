# This script plots a spectrum of an active supermassive black hole
# from Marconi et al. (2004) (https://arxiv.org/pdf/astro-ph/0311619.pdf)
# c.f. Figure 3 for a similar comparison
# the input is the log of the total energy of the active supermassive black hole in ergs^-1


import numpy as np
import matplotlib.pyplot as plt

def Make_SED(log_Lbol_ergs):

    # Convert to Solar luminosity
    log_Lbol = log_Lbol_ergs - 33.58

    # Calculate parameter A_X for normalisation of X-ray part
    A_X = 10**(-13.375 + 1.666*log_Lbol - 0.06283*(log_Lbol)**2 + 0.001402*(log_Lbol)**3)

    # Calculate parameter A_O for normalisation of optical part
    A_O = 10.**((np.log10(A_X) + 2.06)/0.714)

    # Calculate IR and UV normalisations (these aren't free parameters)
    A_IR = A_O*4.74e-36
    A_UV = A_O*2.11e20

    # IR part of spectrum
    log_nu1 = np.linspace(14.25,14.48)
    nuLnu1 = A_IR*(10**log_nu1)**3
    log_nuLnu1 = np.log10(nuLnu1)

    # Optical part of spectrum
    log_nu2 = np.linspace(14.48,15.40,93)
    nuLnu2 = A_O*(10.**log_nu2)**0.56
    log_nuLnu2 = np.log10(nuLnu2)

    # UV part of spectrum
    log_nu3 = np.linspace(15.40,15.78,10)
    nuLnu3 = A_UV*(10.**log_nu3)**(-0.76)
    log_nuLnu3 = np.log10(nuLnu3)

    # X-ray part of spectrum
    log_nu_X = np.array([17.405214248838412, 17.620237480640167, 17.68, 17.781569437274136, 17.956375838926174, 18.131388745482706, 18.225658234383070, 18.279194630872480, 18.387041817243160, 18.535622096024780, 18.630304594734124, 18.711306143520908, 18.819204956117710, 18.980382034073310, 19.141249354672176, 19.288435725348478, 19.435570469798660, 19.542385131646878, 19.622354155911204, 19.728962312854932, 19.808828084667013, 19.875529168817764, 19.928600929272072, 19.981775942178630, 20.048167268972640, 20.101135776974700, 20.154001032524523, 20.220443985544660, 20.273257614868356, 20.339597315436244, 20.37919463087248, 20.445534331440374, 20.498502839442438])

    log_nuLnu_X = np.array([10.188461538461539, 10.207692307692307, 10.21, 10.226923076923077, 10.25, 10.288461538461538, 10.311538461538461, 10.3, 10.334615384615384, 10.403846153846153, 10.457692307692309, 10.492307692307692, 10.53076923076923,  10.538461538461538, 10.523076923076923, 10.488461538461538, 10.45, 10.407692307692308, 10.365384615384615, 10.307692307692307, 10.257692307692308, 10.226923076923077, 10.180769230769231, 10.142307692307693, 10.088461538461537, 10.034615384615385, 9.973076923076922, 9.923076923076923, 9.857692307692307, 9.8,  9.75, 9.692307692307692, 9.638461538461538]) + np.log10(A_X)

    # Interpolating between UV and X-ray
    log_nu4 = np.linspace(15.78, 17.40, 100)
    gradient4 = (log_nuLnu_X[0] - log_nuLnu3[-1])/(17.40-15.78)
    intercept4 = log_nuLnu_X[0] - gradient4*17.40
    log_nuLnu4 = gradient4*log_nu4 + intercept4

    # Put the arrays together
    log_nu = np.concatenate((log_nu1, log_nu2, log_nu3, log_nu4, log_nu_X))
    log_nuLnu = np.concatenate((log_nuLnu1, log_nuLnu2, log_nuLnu3, log_nuLnu4, log_nuLnu_X))
    log_nuLnu = log_nuLnu + 33.58

    # Calculate bolometric luminosity (as a check)
    Lnu = 10.**(log_nuLnu - log_nu)
    Lbol = np.trapz(Lnu, x=(10.**(log_nu)))

    # Calculate alpha OX parameter from 2500 A and 2keV
    alpha_OX_gradient = (log_nuLnu2[60] - log_nu2[60] - log_nuLnu_X[2] + log_nu_X[2])/(17.68 - 15.08)

    # Calculate alpha_OX parameter from 2500 A
    alpha_OX_L = -0.11*(log_nuLnu2[60] - log_nu2[60] + 33.58) + 1.85

    return Lbol, log_nu, log_nuLnu, alpha_OX_gradient, alpha_OX_L


if __name__ == "__main__":
    print 'Making test Spectral Energy Distribution...'
    Lbol1, log_nu1, log_nuLnu1, alpha_OX_gradient1, alpha_OX_L1 = Make_SED(43.)
    Lbol2, log_nu2, log_nuLnu2, alpha_OX_gradient2, alpha_OX_L2 = Make_SED(45.)
    Lbol3, log_nu3, log_nuLnu3, alpha_OX_gradient3, alpha_OX_L3 = Make_SED(47.)

    plt.plot(log_nu1, log_nuLnu1, color='k', ls='solid', label='$L_{\mathrm{bol}} = 10^{43} ergs^{-1}$')
    plt.plot(log_nu2, log_nuLnu2, color='r', ls='dashed', label='$L_{\mathrm{bol}} = 10^{45} ergs^{-1}$')
    plt.plot(log_nu3, log_nuLnu3, color='b', ls='dotted', label='$L_{\mathrm{bol}} = 10^{47} ergs^{-1}$')
    plt.minorticks_on()
    plt.xlabel('$\mathrm{log_{10} (\\nu / Hz)}$', fontsize=14)
    plt.ylabel('$\mathrm{log_{10} (\\nu L_{\\nu} / ergs^{-1})}$', fontsize=14)
    plt.xlim(14.3,21)
    plt.ylim(41,47)
    plt.legend(fontsize='large')

    plt.tick_params(axis='x',which='major',top='on')
    plt.tick_params(axis='x',which='minor',top='on')
    plt.tick_params(axis='y',which='major',right='on')
    plt.tick_params(axis='y',which='minor',right='on')

    plotfile = 'Marconi_SED_test.ps'
    print
    print 'Output plot:'
    print 
    print plotfile
    print 
    plt.savefig(plotfile)





