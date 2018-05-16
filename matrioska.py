import os
import bootmedian as bm
from astropy.io import fits
import numpy as np


def matrioska(target_name, xcen, ycen, incl, PA, nbins, size, c=0,
              mask_name="default_radial.fits", nsimul=1000):

    target_fits = fits.open(target_name)
    radial_mask(target_name=target_name, xcen=xcen, ycen=ycen, incl=incl,
                PA=PA, c=c, mask_name=mask_name)
    mask_fits = fits.open(mask_name)
    r = np.linspace(start=0, stop=size, num=nbins+1)
    r_mids = np.zeros(nbins)
    profile = np.zeros(shape=(8,len(r)), dtype="float32")
    for i in range(len(r_mids)):
        index = [(mask_fits[0].data > r[i]) & (mask_fits[0].data < r[i+1])]
        r_mids[i] = np.mean(mask_fits[0].data[index])
        print(r_mids[i])
        result = [r_mids[i]] + bm.bootmedian(sample_input=target_fits[0].data[index],
                                             nsimul=nsimul, errors=1).tolist()
        profile[:, i] = result

    output = {
              "r_mids": profile[0, :],
              "int_median": profile[1, :],
              "int_s1_up": profile[2, :],
              "int_s1_down": profile[3, :],
              "int_s2_up": profile[4, :],
              "int_s2_down": profile[5, :],
              "int_s3_up": profile[6, :],
              "int_s3_down": profile[7, :]
              }
    return(output)


def radial_mask(target_name, incl, PA, c=0,  xcen=0, ycen=0,
                mask_name="default_radial.fits"):
    target_fits = fits.open(target_name)
    if((xcen == 0) & (ycen == 0)):
        xcen = target_fits[0].header['NAXIS1']/2
        ycen = target_fits[0].header['NAXIS2']/2
    sizex = target_fits[0].header['NAXIS1']
    sizey = target_fits[0].header['NAXIS2']
    q = np.cos(np.radians(incl))
    # def mask_ellipse(mask_name, sizex, sizey, xcen, ycen, q, PA, c):
    mask_ellipse(mask_name=mask_name, sizex=sizex, sizey=sizey, xcen=ycen,
                 ycen=xcen, q=q, PA=PA, c=c)


def mask_ellipse(mask_name, sizex, sizey, xcen, ycen, q, PA, c):
    mask = dist_superellipse((sizex, sizey), (xcen, ycen), q=q, pos_ang=PA, c=c)
    hdu = fits.PrimaryHDU(mask)
    if os.path.exists(mask_name):
        os.remove(mask_name)
    hdulist = fits.HDUList([hdu])
    hdulist[0].data = np.floor(mask)
    hdulist.writeto(mask_name)


def dist_superellipse(n, center, q=1, pos_ang=np.array(0), c=np.array(0)):
    """
    Geometric functions to work with elliptical apertures.

    History:

     24 Aug 2012 - created
     13 Nov 2012 - major restructuring of the teams/miri folder.
              Renamed from mrs_ellipse.py to ellipse.py

     @author: Ruyman Azzollini (CAB, INTA-CSIC)

    Form an array in which the value of each element is equal to the
    semi-major axis of the superellipse of specified center, axial ratio,
    position  angle, and "c" parameter (diskiness/boxiness) which passes
    through that "point".

    Useful for super-elliptical aperture photometry. The term "superellipse",
    a generalization of an ellipse that may be "disky" (c<0), or "boxy" (c>0),
    is borrowed from Peng et al. 2002 (GALFIT).

    Inspired on dist_ellipse.pro from AstroLib (IDL), but

    WARNING: this program doesn't take into account the change in the order of
    axes from IDL to Python. That means, that in 'n' and in 'center', the order
    of the coordinates must be reversed with respect to the case for
    dist_ellipse.pro, in order to get expected results. Nonetheless, the polar
    angle means the counter-clock wise angle with respect to the 'y' axis.

    EXPLAIN BETTER!!
    Finish Documentation
    """

    print(q)
    print(pos_ang)
    print(c)

    criterium1 = type(n) in (type(1), type(1.0)) or \
    type(n) == type(()) and len(n) == 2

    if not criterium1:
        print('n must be an scalar or a 2 elements tuple')

    if type(n) in (type(1), type(1.0)):
        n = float(n)
    if type(n) == 'tuple':
        n = tuple([float(x) for x in n])

    criterium2 = 'int' in str(type(q)) or 'float' in str(type(q))
    if not criterium2:
        print 'q must be an integer or float'
        criterium3 = False
    else:
        criterium3 = q >= 0 and q <= 1

    if not criterium3:
        print('q = %f is out of range (0,1)' % q)

    if criterium2 and criterium3:
        q = float(q)

    criterium4 = 'int' in str(type(pos_ang)) or 'float' in str(type(pos_ang))
    if not criterium4:
        print 'pos_ang must be an integer or float'
        criterium5 = 'False'
    else:
        criterium5 = -180 <= pos_ang <= 180

    if not criterium5:
        print 'pos_ang out of range (-180,180)'

    if criterium4 and criterium5:
        pos_ang = float(pos_ang)

    criterium6 = 'int' in str(type(c)) or 'float' in str(type(c))
    if not criterium6:
        print 'c must be an integer or float'
    else:
        c = float(c)

    criterium7 = type(center) == type(()) and len(center) == 2
    if not criterium7:
        print 'center must be a 2 element tuple'
    else:
        center = tuple([float(x) for x in center])

    # END CHECK OF INPUTS

    criteriums = reduce(np.logical_and,np.array([criterium1,criterium2,criterium3,\
    criterium4,criterium5,criterium6,criterium7]))

    if criteriums :

        radeg = 180. / np.pi

        ang = pos_ang / radeg
        cosang = np.cos(ang)
        sinang = np.sin(ang)

        if len(n) == 2 :
            nx = n[1]
            ny = n[0]
        else :
            nx = ny = n

        x = np.arange(nx,dtype='Float32') - center[1]
        y = np.arange(ny,dtype='Float32') - center[0]
        im = np.zeros(shape=(ny,nx),dtype='Float32')
        xcosang = x * cosang
        xsinang = x * sinang

        for i in range(ny) :
            xtemp = xcosang + y[i] * sinang
            ytemp = -xsinang + y[i] * cosang
            im[i,:] = ( (np.abs(xtemp/q))**(c+2.) + \
			(np.abs(ytemp))**(c+2.) )**(1./(c+2.))

    else: raise RuntimeError

    return im
