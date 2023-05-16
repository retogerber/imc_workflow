# Usage

1. Move .mcd files to directory results/raw

2. Update file config/channel_metadata.csv which contains three columns:
    - channel: Identifier, like isotope plus mass
    - name: Name of channel, used for matching
    - type: one of 'type' (cell type marker), 'state' (cell state marker), 'seg' (segmentation marker or 'unknown' (other channels)

3. Check/Update config/config.yaml
    - filtering:
        - `filter_thres`
    - segmentation:
        - `nuc_channels`: names (matching column `name` in config/channel_metadata.csv ) of nuclear markers used for segmentation
        - `mem_channels`: names (matching column `name` in config/channel_metadata.csv ) of membrane markers used for segmentation

4. Run the workflow with 
```
snakemake --use-conda --use-singularity -c 1
```
the flag `-c` allows you to set the number of cores to use

5. (Optional) Visualize segmentation
First create environment:
```
conda env create -f workflow/env/napari_vis.yaml
```
then run (replace IMAGE_NAME with a name from config/samples.csv)
```
conda run -n napari_vis python workflow/scripts/visualize_mask.py -s IMAGE_NAME
```
to open napari and visually inspect the segmentation. To show the filtered image use the flag `-f` in the above command.

If the segmentation is insufficient you can either change the filtering or change mesmer postprocessing parameters (file config/mesmer_postprocess_id01.yml, see [here](https://deepcell.readthedocs.io/en/master/_modules/deepcell/applications/mesmer.html)).

6. Update file config/sample_metadata.csv:
    It should contain at least two columns, the first one has to be `sample` the image names (check the file config/samples.csv for all available names). Other columns are used to describe the individual samples.

7. Update `cols_for_batch` parameter in config/config.yaml, add all batch relevant column names from config/sample_metadata.csv. All columns used for batch correction should be categorical. The first one should be `sample`.

8. Run the workflow with 
```
snakemake --use-conda --use-singularity -c 1 
```

9. Create and inspect report
To create a report run:
```
snakemake --use-conda --use-singularity -c 1 --report
```
then open the file report.html . In the top left are links to reports from the individual steps (`cell level filtering report`, `batch correction report`, ` batch correction evaluation report`, `cell clustering report`). If filtering and batch correction is insufficient either change parameters in `config/config.yaml` or switch to more manual analysis.

10. Annotation
Choose a clustering in the report from above (link `cell clustering report`) and manually annotate the clusters. Add the annotation to file config/clustering_annotation.csv which is a two column csv file. The first column has name `cluster` and contains the numbers of the choosen cluster. The second column has name `celltype` and contains the annotations. Additionally in config/config.yaml update parameter `clustering_used` to the cluster name that has been used for annotation.

11. add annotation
Run again
```
snakemake --use-conda --use-singularity -c 1 
```
The final output file is results/spe/filt/SPE_combined_03.rds

