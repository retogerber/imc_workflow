import numpy as np
from skimage.io import imread
from skimage.io import imsave
from scipy.ndimage import maximum_filter

img_name_in = snakemake.input[0]
img_name_out = snakemake.output[0]
thres = snakemake.params[0]
subtraction_type = snakemake.params[1]
#thres = 1

def background_subtraction(img: np.ndarray, thres: float) -> np.ndarray:
    #img = img.astype(int)
    if subtraction_type == "below_threshold":
        img = np.where(img <= thres, 0, img)
    elif subtraction_type == "all":
        img = np.where(img <= thres, 0, img-thres)
    else:
        raise(Exception("invalid subtraction_type"))
    #return img.astype(np.uint16)
    return img

# read image
img = imread(img_name_in)

# apply hot pixel filter
img_filt = background_subtraction(img, thres)

# save image
imsave(img_name_out, img_filt)
