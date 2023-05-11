# General configuration

1. place the .mcd files in results/raw
2. update the file config/channel_metadata.csv
    the file contains three columns: channel, name, type
        channel: channel ID (e.g. element + mass)
        name: name of the channel, used for matching
        type: one of 
            "seg", segmentation marker, used for segmentation with mesmer
            "type", cell type marker, used for clustering and cell type annotation
            "state", cell state marker, other markers not used in this workflow
        
3. update the main config file config/config.yaml
    important points:
        specify the nuclear channels, has to match entry from "name" column from file config/channel_metadata.csv (parameter "segmentation.nuc_channels")
        specify the membrane channels, has to match entry from "name" column from file config/channel_metadata.csv (parameter "segmentation.mem_channels")
        specify the filtering type (parameter "filtering.filter_subtraction_type") and the threshold (parameter "filtering.filter_thres")
4. (optional) change postprocessing steps of mesmer:
    change file config/mesmer_postprocess_id01.yml
    for options see [here](https://deepcell.readthedocs.io/en/master/_modules/deepcell/applications/mesmer.html)

5. Run workflow with
```
snakemake --use-conda --use-singularity
```

6. (optional) Check visually the quality of the segmentation
    First create a virtual environment using conda (or mamba):
```
mamba env create -f workflow/envs/napari_vis.yaml
```
Then visualize the created mask overlayed on the image using napari by running (replace SAMPLE_NAME with a valid sample name from the file config/samples.csv):
```
conda run -n napari_vis python workflow/scripts/visualize_mask.py -s SAMPLE_NAME
```
To show the filtered image (if you applied a filtering step) add the "-f" flag to the above command.


7. Update file config/sample_metadata.csv
   used for batch correction, needs at least two columns. The first column should be named "sample" and should contain the sample names (see the newly created file config/samples.csv). Additional columns should contain information about the individual samples.
    The batch effect to remove should be specified with the parameter "cols_for_batch" in the config/config.yaml file

8. Rerun workflow with
```
snakemake --use-conda --use-singularity
```
9. (optional) Check reports 
    check the output reports in the directory results/html

10. Manual annotation
    using html report results/html/visualize_cluster_spe.html check the clustering result. Choose a clustering which best represents celltypes and add the clustering name to the parameter "clustering_used" in the config/config.yaml file. Then write down the manual annotation in the file config/clustering_annotation_file.csv which consists of two columns: cluster (the cluster number) and celltype (the name of the cluster) 

11. Rerun workflow with
```
snakemake --use-conda --use-singularity
```

The final output file is results/spe/filt/SPE_combined_03.rds which is an R object of class SpatialExperiment that can be used for further downstream analysis.


