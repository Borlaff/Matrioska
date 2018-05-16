# Matrioska:
## A program to create surface brightness profiles in Python
### Author: Alejandro Serrano Borlaff
### Instituto de Astrofisica de Canarias
### asborlaff@gmail.com

### v1.0 - 16 May 2018 - First working version.

If you use this software, please cite Borlaff et al. 2018:

"Evolution of the anti-truncated stellar profiles of S0 galaxies since z=0.6 in the SHARDS survey: II - Structural and photometric evolution" - Astronomy & Astrophysics



 Input:
 
    target_name: Fits to be analysed
    
    xcen and ycen: Coordinates of the center of the object.
    
    incl: Inclination of the object (incl = np.cos^-1(b/a))
    
    PA: Position angle of the mayor axis, measured from the positive y-axis.
    
    size: Max radius to calculate the profile.
    
    c: Generatization of ellipse, "disky" (c<0), or "boxy" (c>0): Default=0.
    
    mask_name : Name of the fits file that contain the elliptical apertures.
    
    nsimul: Number of bootstrapping simulations for bootmedian
    
 Output: Python Dictonary:
 
               "r_mids": Central pixel position of the elliptical aperture. 
               
              "int_median": Median intensity of the elliptical aperture
              
              "int_s1_up": 1 sigma upper limit of the probability distribution of median intensity.
              
              "int_s1_down": 1 sigma upper limit of the probability distribution of median intensity.
              
              "int_s2_up": Same for 2 sigma level. 
              
              "int_s2_down": Same for 2 sigma level. 
              
              "int_s3_up": Same for 3 sigma level.
              
              "int_s3_down": Same for 3 sigma level. 

