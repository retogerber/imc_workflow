import numpy as np
from skimage.io import imread
from skimage.io import imsave
from skimage.filters import threshold_minimum
from scipy.ndimage import distance_transform_edt
from cv2 import filter2D
import diptest

img_name_in = snakemake.input[0]
img_name_out = snakemake.output[0]
thres = snakemake.params["thres"]
subtraction_type = snakemake.params["subtraction_type"]
kernelsize = snakemake.params["kernelsize"]

def convolve_filter_image(img, kersize):
    # create kernel
    kersize = kersize-1 if kersize%2==0 else kersize
    kernel = np.ones((kersize,kersize))
    kernel[int((kersize-1)/2),int((kersize-1)/2)] = 0
    # distance from center pixel
    distkernel = distance_transform_edt(kernel)
    distkernel[int((kersize-1)/2),int((kersize-1)/2)] = 1
    # invert distance
    kernel = 1/distkernel
    # set max distance, i.e. disk kernel
    kernel[kernel<kernel[int((kersize-1)/2),0]]=0
    kernel[int((kersize-1)/2),int((kersize-1)/2)] = 1
    # run convolution
    tmp = filter2D(src = img, ddepth=-1, kernel=kernel)

    # test if bimodal
    _, pval = diptest.diptest(tmp.flatten())

    # not bimodal
    if pval>0.05:
        thresh = np.quantile(tmp, 0.01)
        img[tmp<thresh] = 0

    # bimodal
    else:
        otsuin = tmp[img>0]
        otsuin = otsuin[otsuin<np.quantile(otsuin, 0.9)]
        try:
            thresh = threshold_minimum(otsuin)

        # if cannot find two peaks
        except RuntimeError:
            thresh = np.quantile(otsuin, 0.01)
        nmax=1
        n=0
        npix=[np.sum(tmp>thresh)+1,np.sum(tmp>thresh)]
        while n<nmax and npix[-1]!=npix[-2]:
            n+=1
            img[tmp<thresh] = 0
            tmp = filter2D(src = img, ddepth=-1, kernel=kernel)
            npix.append(np.sum(tmp>thresh))
    return img

def background_subtraction(img: np.ndarray, thres: float, kernelsize: float) -> np.ndarray:
    #img = img.astype(int)
    if subtraction_type == "below_threshold":
        img = np.where(img <= thres, 0, img)
    elif subtraction_type == "all":
        img = np.where(img <= thres, 0, img-thres)
    elif subtraction_type == "convolve":
        for i in range(img.shape[0]):
            img[i,:,:] = convolve_filter_image(img[i,:,:].copy(), kernelsize) 
    else:
        raise(Exception("invalid subtraction_type"))
    #return img.astype(np.uint16)
    return img

# read image
img = imread(img_name_in)

# apply hot pixel filter
img_filt = background_subtraction(img, thres, kernelsize)

# save image
imsave(img_name_out, img_filt)
