import sys, os
import sav2fits as sf
import bootmedian as bm
import ellipse as ell
from astropy.io import fits
import numpy as np
from tqdm import *
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel, Gaussian2DKernel
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
import matplotlib.pyplot as plt
# Define model
@custom_model


def gaussian_polar(x, amplitude=1., mean=1., sd=1.):    
    order = 10
    x_pred = np.linspace(start=-order*360,stop=(order+1)*360,num=(2*order+1)*360+1)
    while mean<0:
        mean = mean + 360
    while mean>360:
        mean = mean - 360
    gauss1D = models.Gaussian1D(amplitude=amplitude,mean=mean,stddev=sd)
    y_pred =  gauss1D(x_pred)
    index_cen = ((x_pred > 0) & (x_pred < 360)) 
    #x_pred = x_pred[index_cen]
    y_out = y_pred[index_cen]
    
    for i in np.linspace(start=1,stop=order, num=order):
        index_up = ((x_pred > i*360) & (x_pred < (i+1)*360))
        index_down = ((x_pred < -(i-1)*360) & (x_pred > -i*360)) 
        y_out = y_out + y_pred[index_up] + y_pred[index_down]
    y = np.interp(x=x, xp=x_pred[index_cen], fp=y_out)
    return(y)



def find_spiral(image_name, PA, incl, xcen=0, ycen=0, rmin=10, rmax=0):
    m0 = fits.open(image_name)

    kernel = Gaussian2DKernel(stddev=1)
    intensity = convolve(m0[0].data, kernel)
    #intensity = m0[0].data
    if (xcen == 0):
        xcen = m0[0].header['NAXIS1']/2
    if (ycen == 0):
        ycen = m0[0].header['NAXIS2']/2
    if (rmax == 0):
        rmax = xcen
    
    angle_mask(target_name=image_name, c=0,incl=incl,PA=PA, mask_name=image_name.replace(".fits","")+"angle_mask.fits", xcen = xcen, ycen=ycen)
    radial_mask(target_name=image_name, c=0,incl=incl,PA=PA, mask_name=image_name.replace(".fits","")+"radial_mask.fits", xcen = xcen, ycen=ycen)

    radial = fits.open(image_name.replace(".fits","") + "radial_mask.fits")
    angle = fits.open(image_name.replace(".fits","") + "angle_mask.fits")
    
    r = radial[0].data
    theta = angle[0].data

    print(rmax)

    r_bins = np.linspace(start=rmin,stop=rmax,num=100)
    
    i = range(len(r_bins)-1)
    
    radii = np.zeros(len(r_bins)-1)
    mean_0 = np.zeros(len(r_bins)-1)
    sd_0 = np.zeros(len(r_bins)-1)
    amp_0 = np.zeros(len(r_bins)-1)
    mean_1 = np.zeros(len(r_bins)-1)
    sd_1 = np.zeros(len(r_bins)-1)
    amp_1 = np.zeros(len(r_bins)-1)
                    
    for i in tqdm(i):
        index = ((radial[0].data > r_bins[i]) & (radial[0].data < r_bins[i+1]))
        y = intensity[index]
        x = angle[0].data[index]
        xy = zip(*sorted(zip(x,y)))
        x = np.array(xy[0])
        y = np.array(xy[1])
        bool_isnan = (np.isnan(y))
        x = x[~bool_isnan]
        y = y[~bool_isnan]
        y = y - np.nanpercentile(y,5)

        amplitude_1 = np.max(y[(x < 180)])
        amplitude_2 = np.max(y[(x > 180)])

        if (i == 0):
            bool_peak_1 = (y == amplitude_1)
            bool_peak_2 = (y == amplitude_2)
            peak_1 = x[bool_peak_1][0]
            peak_2 = x[bool_peak_2][0]
            sd = 20.

        if (i > 0):
            peak_1 = mean_0[i-1]
            peak_2 = mean_1[i-1]
#            amplitude_1 = amp_0[i-1]
#            amplitude_2 = amp_1[i-1]
            sd = sd_0[i-1]
        print(peak_1)
        print(peak_2)
        
        model_2gauss_polar = gaussian_polar(amplitude=amplitude_1, mean=peak_1, sd=sd) + gaussian_polar(amplitude=amplitude_2, mean=peak_2, sd=sd)
        fitter = fitting.SLSQPLSQFitter()
        gg_fit = fitter(model_2gauss_polar, x, y)
        radii[i] = np.median(r[index])
        mean_0[i] = gg_fit.mean_0[0]
        mean_1[i] = gg_fit.mean_1[0]
        amp_0[i] = gg_fit.amplitude_0[0]
        amp_1[i] = gg_fit.amplitude_1[0]
        sd_0[i] = gg_fit.sd_0[0]
        sd_1[i] = gg_fit.sd_1[0]

        plt.figure(figsize=(8,5))
        plt.plot(x, y, 'ko')
        plt.plot(x, gg_fit(x))
        plt.xlabel('Position')
        plt.ylabel('Flux')
        plt.show()   

    return([radii, mean_0, mean_1, amp_0, amp_1, sd_0, sd_1])

def axis_vel(m1_name, PA, incl, vsys, rmax, xcen=0, ycen=0):
    
    #path = "/home/borlaff/PHD/NGC3433/"+name_sim+"/"
    # We create the fake obsertavions from the simulation 
    #sf.sav2cube(path=path, time=time, type_p = type_p, xunit = xunit, vunit=vunit, pa = PA, incl = incl, 
    #            pixscale=pixscale, transpose=True, dist=dist, sd_numeric=sd_numeric, vsys=vsys)
    #rotcur = sf.get_rotcur(dist=dist, path=path, time=time, type_p=type_p, vunit=vunit, xunit=xunit)

    #m0 = name_sim+type_p+"xyz"+str(time)+"_m0_PA"+str(np.round(np.degrees(np.radians(PA)),0))+"_i"+str(np.round(np.degrees(np.radians(incl)),0))+".fits"
    #m1 = name_sim+type_p+"xyz"+str(time)+"_m1_PA"+str(np.round(np.degrees(np.radians(PA)),0))+"_i"+str(np.round(np.degrees(np.radians(incl)),0))+".fits"
    # We open the fits with the fake observations 
    #m0_fits = fits.open(m0_name) # Intensity 2D field
    m1_fits = fits.open(m1_name) # Velocity 2D field
    if((xcen == 0) & (ycen==0)): 
        xcen = m1_fits[0].header['NAXIS1']/2
        ycen = m1_fits[0].header['NAXIS2']/2

    angle_mask(target_name=m1_name,xcen=xcen, ycen=ycen, PA=PA,incl=incl,c=0,mask_name=m1_name.replace(".fits","")+"_angle_mask.fits")
    radial_mask(target_name=m1_name,xcen=xcen, ycen=ycen, PA=PA,incl=incl,c=0,mask_name=m1_name.replace(".fits","")+"_radial_mask.fits")
    # profile=matrioska.matrioska(target_name="ngc3433_839sxyz"+str(time)+"_m0_PA22.0_i37.0.fits", xcen=xcen,ycen=ycen, 
    #                             incl = incl, PA = PA, nbins=50, size =  xcen, c=0)
    angle_fits = fits.open(m1_name.replace(".fits","")+"_angle_mask.fits") # Radial mask
    radial_fits = fits.open(m1_name.replace(".fits","")+"_radial_mask.fits") # Angle mask

    # We select the pixels in the PA line and the Node line
    wide = 10
    PA_index = ((np.abs(angle_fits[0].data - PA) < wide) | (np.abs(angle_fits[0].data - PA - 180) < wide))
    Node_index = ((np.abs(angle_fits[0].data - PA + 90) < wide) | (np.abs(angle_fits[0].data - PA - 90) < wide))
    radial_index = (radial_fits[0].data < rmax)
    

    r_PA = radial_fits[0].data[(PA_index & radial_index)] 
    v_PA = np.abs(m1_fits[0].data[(PA_index & radial_index)] - vsys)
    r_Node = radial_fits[0].data[(Node_index & radial_index)]
    v_Node = m1_fits[0].data[(Node_index & radial_index)] - vsys
    
    r_node_median, v_node_median = median_line(x=r_Node, y=v_Node, nbins=20)
    r_PA_median, v_PA_median = median_line(x=r_PA, y=v_PA, nbins=20)
    
    return([[r_node_median, v_node_median],[r_PA_median, v_PA_median]])

def median_line(x,y,nbins):

    x_bins = np.linspace(start=np.nanmin(x), stop=np.nanmax(x), num=nbins+1)
    x_mids    = np.zeros(nbins)
    x_1s_up   = np.zeros(nbins)
    x_1s_down = np.zeros(nbins)
    y_mids    = np.zeros(nbins)
    y_1s_up   = np.zeros(nbins)
    y_1s_down = np.zeros(nbins)
    for i in range(nbins):
        try:
            index = ((x >= x_bins[i]) & (x < x_bins[i+1]))
            y_bm = bm.bootmedian(y[index], nsimul=10000, errors=1)
            x_bm = bm.bootmedian(x[index], nsimul=10000, errors=1)

        except:
            print("All nan")
            y_bm = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
            x_bm = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
        x_mids[i] = x_bm[0]
        x_1s_up[i] = x_bm[1]
        x_1s_down[i] = x_bm[2]

        y_mids[i] = y_bm[0]
        y_1s_up[i] = y_bm[1]
        y_1s_down[i] = y_bm[2]
    return([x_mids,x_1s_up, x_1s_down, y_mids, y_1s_up, y_1s_down])
    

def matrioska(target_name, xcen, ycen, incl, PA, nbins, size, c=0, mask_name = "default_radial.fits"):
    target_fits = fits.open(target_name)
    sizex = target_fits[0].header['NAXIS1']
    sizey = target_fits[0].header['NAXIS2']
    radial_mask(target_name, xcen, ycen, incl, PA, c, mask_name=mask_name)
    angle_mask(target_name, xcen, ycen, incl, PA, c)
    mask_fits = fits.open(mask_name)
    #max_rad = np.max(mask_fits[0].data)
    r = np.linspace(start=0, stop=size, num=nbins+1)
    r_mids = np.zeros(nbins)
    profile = np.zeros(shape=(8,len(r)), dtype="float32")
    for i in range(len(r_mids)):
        index = [(mask_fits[0].data > r[i]) & (mask_fits[0].data < r[i+1])]
        r_mids[i] = np.mean(mask_fits[0].data[index])
        print(r_mids[i])
        #print([r[i]] + bm.bootmedian(sample_input=target_fits[0].data[index], nsimul=10000, errors=1))
        result = [r_mids[i]] + bm.bootmedian(sample_input=target_fits[0].data[index], nsimul=10000, errors=1).tolist()
        profile[:,i] = result
    return(profile)




def radial_mask(target_name, incl, PA, c=0,  xcen=0, ycen=0, mask_name = "default_radial.fits"):
    target_fits = fits.open(target_name)
    if((xcen == 0) & (ycen==0)): 
        xcen = target_fits[0].header['NAXIS1']/2
        ycen = target_fits[0].header['NAXIS2']/2    
    sizex = target_fits[0].header['NAXIS1']
    sizey = target_fits[0].header['NAXIS2']
    q = np.cos(np.radians(incl))
    # def mask_ellipse(mask_name, sizex, sizey, xcen, ycen, q, PA, c):

    a=ell.mask_ellipse(mask_name = mask_name,sizex = sizex, sizey = sizey, xcen=ycen, ycen=xcen, q=q, PA=PA, c=c)

def angle_mask(mask_name, incl, PA, xsize, ysize, c=0, xcen=0, ycen=0):
    if((xcen == 0) & (ycen==0)): 
        xcen = xsize/2
        ycen = ysize/2 

    q = np.cos(np.radians(incl))
    a=ell.mask_ellipse(mask_name, ysize, xsize, ycen, xcen, q, PA, c)
    mask_fits = fits.open(mask_name)
    # r^2 = x^2 + y^2
    Ax_array = np.linspace(start=0, stop=xsize-1, num=xsize) - xcen
    Ay_array = np.linspace(start=0, stop=ysize-1, num=ysize) - ycen
    
    for i in range(xsize):
    	for j in range(ysize):
            mask_fits[0].data[j,i] = np.degrees(np.arctan2(Ay_array[j],Ax_array[i]) - np.pi/2)
    while (((mask_fits[0].data<0) - (mask_fits[0].data>360)).any()):
        below_zero = (mask_fits[0].data<0)
        over_360 = (mask_fits[0].data>360)
        if((below_zero).any()): mask_fits[0].data[below_zero] = mask_fits[0].data[below_zero] + 360
        if((over_360).any()): mask_fits[0].data[over_360] = mask_fits[0].data[over_360] - 360

    mask_fits.verify("silentfix")
    os.system("rm "+mask_name)
    mask_fits.writeto(mask_name)


def get_nodeline(vfield, xcen, ycen, PA, incl, size, nbins=30, vsys=0):
    radial_mask(target_name=vfield, xcen=xcen, ycen=ycen, incl=incl, PA=PA, c=0)
    angle_mask(target_name=vfield, xcen=xcen, ycen=ycen, incl=incl, PA=PA, c=0)
    radial_fits = fits.open("default_radial.fits") # Radial mask
    angle_fits = fits.open("default_angle.fits") # Angle mask 
    m1_fits = fits.open(vfield) # Velocity 2D field 

    bins = np.linspace(start=0, stop=size, num=nbins) # We create a mask 
    mids = bins[0:-1]
    westline = []
    westline_up = []
    westline_down = []    
    eastline = []
    eastline_up = []
    eastline_down = []    
    for i in range(len(bins)-1):
        index = ((radial_fits[0].data > bins[i]) & (radial_fits[0].data < bins[i+1]))
        x = angle_fits[0].data[index]
        y = m1_fits[0].data[index] - vsys

        bool_isnan = np.isnan(y)
        x = x[~bool_isnan]
        y = y[~bool_isnan]

        nodeline = get_node(x,y,PA)
        westline.append(nodeline[0][0])
        westline_up.append(nodeline[0][1])
        westline_down.append(nodeline[0][2])

        eastline.append(nodeline[1][0])
        eastline_up.append(nodeline[1][1])
        eastline_down.append(nodeline[1][2])

    return([np.array(mids),np.array(bins),np.array(westline),np.array(westline_up), np.array(westline_down),
     np.array(eastline), np.array(eastline_up), np.array(eastline_down)])
# Get_node uses the distance to the PA+90 and PA+270
# to get the median angle of the line on nodes, at a certain radius. 
def get_node(x, y, PA):
    distance_to_PA = np.degrees(np.arctan2(np.sin(np.radians(x-PA)), np.cos(np.radians(x-PA))))
    distance_to_PA[distance_to_PA < 0] = distance_to_PA[distance_to_PA < 0] + 360
    west_index = (distance_to_PA > 0) & (distance_to_PA < 180) 
    east_index = (distance_to_PA > 180) & (distance_to_PA < 360)
    #print(distance_to_PA[west_index])
    #print(distance_to_PA[east_index])
#    print(x[west_index])
#    print((1/y[west_index]**2))
    node_west = bm.bootmedian(sample_input = x[west_index], nsimul=10000, errors=1, weights_input = (1./(abs(y[west_index]))))
    node_east = bm.bootmedian(sample_input = x[east_index], nsimul=10000, errors=1, weights_input = (1./(abs(y[east_index]))))
    #plt.plot(x, y, 'ro')
    #plt.plot(node_west, 0, "bo")
    #plt.plot(node_east, 0, "bo")
    #plt.show()
    return([node_west, node_east])
