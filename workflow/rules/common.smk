import numpy as np

def get_samples_from_images_csv(wildcards):
    """Return samples from file 'results/images.csv'"""
    with checkpoints.create_image_csv.get(**wildcards).output[0].open() as f:
        samples_df = pd.read_csv(f)
        images = samples_df["image"].tolist()
        samples = [ os.path.splitext(im)[0] for im in images ]
        return samples

def get_img_from_images_csv(wildcards):
    samples = get_samples_from_images_csv(wildcards)
    samples = np.array(samples)[np.array(samples)==wildcards.sample]
    samples_ls = [ f"results/img/{sam}.tiff" for sam in samples ]
    return samples_ls

def get_images_from_images_csv(wildcards):
    samples = get_samples_from_images_csv(wildcards)
    samples_ls = [ f"results/img/{sam}.tiff" for sam in samples ]
    return samples_ls

def mcd_name_from_sample_name(wildcards):
    str_output = re.sub("_[0-9]{3}$", "", wildcards.sample)
    mcd_name = f"results/raw/{str_output}.mcd"
    return mcd_name

def mcd_summary_panel_from_sample_name(wildcards):
    str_output = re.sub("_[0-9]{3}$", "", wildcards.sample)
    mcd_summary_panel = f"results/summary_panels/{str_output}_summary.csv"
    return mcd_summary_panel


def get_img_samplenames(wildcards):
    '''return full path tiff from wildcard'''
    out = []
    for s in wildcards:
        if config["filtering"]["filter_thres"] > 0:
            out.append(f"results/img_filt/{s}.tiff")
        else:
            out.append(f"results/img/{s}.tiff")
    return out

def get_spe_from_images_csv(wildcards):
    samples = get_samples_from_images_csv(wildcards)
    samples_ls = [ f"results/spe/raw/SPE_raw_{sam}.rds" for sam in samples ]
    return samples_ls


def get_mask_from_images_csv(wildcards):
    samples = get_samples_from_images_csv(wildcards)
    samples_ls = [ f"results/masks/{sam}.tiff" for sam in samples ]
    return samples_ls


