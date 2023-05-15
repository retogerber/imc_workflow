
def get_samples_from_images_csv(wildcards):
    """Return samples from file 'results/images.csv'"""
    with checkpoints.create_image_csv.get(**wildcards).output[0].open() as f:
        samples_df = pd.read_csv(f)
        images = samples_df["image"].tolist()
        samples = [ os.path.splitext(im)[0] for im in images ]
        return samples

def get_img_from_images_csv(wildcards):
    samples = get_samples_from_images_csv(wildcards)
    return [ f"results/img/{sam}.tiff" for sam in samples ]

#rule create_samples_csv:
#    input:
#        imcsv="results/images.csv"
#    output:
#        samcsv="config/samples.csv"
#    shell:
#        "cat {input.imcsv} | cut -d, -f1 | sed 1d | sed 's/.tiff//g' > {output.samcsv}"

def mcd_name_from_sample_name(wildcards):
    str_output = re.sub("_[0-9]{3}$", "", wildcards.sample)
    return f"results/raw/{str_output}.mcd"

def mcd_summary_panel_from_sample_name(wildcards):
    str_output = re.sub("_[0-9]{3}$", "", wildcards.sample)
    return f"results/summary_panels/{str_output}_summary.csv"


def get_img_samplenames(wildcards):
    '''return full path tiff from wildcard'''
    out = []
    for s in wildcards:
        if filter_thres > 0:
            out.append(f"results/img_filt/{s}.tiff")
        else:
            out.append(f"results/img/{s}.tiff")
    return out

def get_spe_from_images_csv(wildcards):
    samples = get_samples_from_images_csv(wildcards)
    return [ f"results/spe/raw/SPE_raw_{sam}.rds" for sam in samples ]


def get_mask_from_images_csv(wildcards):
    samples = get_samples_from_images_csv(wildcards)
    return [ f"results/masks/{sam}.tiff" for sam in samples ]


