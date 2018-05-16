import numpy as np
import multiprocessing
import bottleneck as bn

def median_bootstrap(argument):
    # arguments = sample, indexes, i
    sample = argument[0]
    try:
    	weigths = argument[1]
    except: 
    	weigths = (np.zeros(len(sample))+1)
    weigths = weigths/np.sum(weigths)
    x = np.random.choice(sample,size=len(sample),p=weigths)   
    median_boot = bn.nanmedian(x)
    return(median_boot)

def bootmedian(sample_input, nsimul, errors=1, weights_input = False):
    #print("Bootmedian v.3.0")
    sigma1 = 0.682689492137086
    sigma2 = 0.954499736103642
    sigma3 = 0.997300203936740

    if(len(sample_input) == 0):
        return(np.array([0,0,0,0,0,0,0]))
    if(len(sample_input) == 1):
        return(np.array([sample_input[0],0,0,0,0,0,0]))


    s1_down_q = (1-sigma1)/2 
    s1_up_q = 1 - s1_down_q
    s2_down_q = (1-sigma2)/2 
    s2_up_q = 1 - s2_down_q
    s3_down_q = (1-sigma3)/2 
    s3_up_q = 1 - s3_down_q
   
    sample_0 = np.ndarray.flatten(sample_input)
    sample = sample_0[~(np.isnan(sample_0))]
    zip_sample = [sample]*nsimul
    median_boot = np.empty(nsimul)
    num_cores = multiprocessing.cpu_count() - 2
    n = len(sample)
    # print("A total of "+str(num_cores)+" workers joined the cluster!")

    if not isinstance(weights_input, bool): 
        weigths_0 = np.ndarray.flatten(weights_input)
        weigths = weigths_0[~(np.isnan(weigths_0))]
    	zip_weigths = [weigths]*nsimul
        arguments = zip(zip_sample,zip_weigths)
    else:
        arguments = zip(zip_sample)


    #return(arguments)
    pool = multiprocessing.Pool(processes=num_cores)
    median_boot = pool.map(median_bootstrap, arguments)
    #print(median_boot)
    pool.terminate()

    median = bn.nanmedian(median_boot)
    if(errors == 1):
        s1_up = np.percentile(median_boot,s1_up_q*100)
        s1_down = np.percentile(median_boot,s1_down_q*100)
        s2_up = np.percentile(median_boot,s2_up_q*100)
        s2_down = np.percentile(median_boot,s2_down_q*100)
        s3_up = np.percentile(median_boot,s3_up_q*100)
        s3_down = np.percentile(median_boot,s3_down_q*100)

    if(errors == 0):
        s1_up = 0
        s1_down = 0
        s2_up = 0
        s2_down = 0
        s3_up = 0
        s3_down = 0

    return(np.array([median, s1_up, s1_down, s2_up, s2_down, s3_up, s3_down]))
