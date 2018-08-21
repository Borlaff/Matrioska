# Demo for matrioska
# ------------------------------------#
# Alejandro Serrano Borlaff
# Instituto de Astrofisica de Canarias
# asborlaff@gmail.com
#######################################


import matrioska as ma

profile = ma.matrioska(target_name="MESSIER_101_I_103aE_dss1_cropped.fits", target_ext=1,
                       xcen=160, ycen=111, incl=60, PA=90, nbins=100,
                       size=150, c=0)
