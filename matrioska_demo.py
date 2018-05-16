# Demo for matrioska
# ------------------------------------#
# Alejandro Serrano Borlaff
# Instituto de Astrofisica de Canarias
# asborlaff@gmail.com
#######################################


import matrioska as ma

profile = ma.matrioska(target_name="MESSIER_101_I_103aE_dss1.fits",
                       xcen=458, ycen=458, incl=30, PA=45, nbins=100,
                       size=150, c=0)
