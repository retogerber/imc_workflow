from tifffile import imread
import cv2
import napari
import pandas as pd
import yaml
import numpy as np
import click

if False:
    import os
    os.chdir("Nextcloud/Projects/imc_workflow")


@click.command(name="visualize_mask", help="open napari to visualize mask")
@click.option("-s","--sample")
@click.option("-i","--img", default=None)
@click.option("-m","--mask", default=None)
@click.option("-p","--panel", default=None)
@click.option("-c","--config",default="config/config.yaml")



def visualize(sample, img, mask, panel, config):
    """Reads in sample and opens napari for visualization """ 
    
    if img is None:
        img = "results/img/"+sample+".tiff"
    if mask is None:
        mask = "results/masks/"+sample+".tiff"
    if panel is None:
        panel = "results/summary_panels/"+sample+"_summary.csv"    
    
    # load summary panel, extract channel names
    sumdf = pd.read_csv(panel)
    channels = sumdf["name"].tolist()
    # load config file to get nuclear and membrane channels
    with open(config, 'r') as stream:
        try:
            parsed_yaml=yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            
    # load image and reorder channels
    im=imread(img)
    
    # channel number of nuclear markers
    nucs=[]
    for c in parsed_yaml["segmentation"]["nuc_channels"]:
        nucs.append(np.where([x==c for x in channels])[0][0])
    
    # channel number of membrane markers
    seg=[]
    for c in parsed_yaml["segmentation"]["mem_channels"]:
        seg.append(np.where([x==c for x in channels])[0][0])
    
    
    # create new image layer order so that nuclear and membrane channels are at the top
    c_ord = seg + nucs
    for i in range(len(channels)):
        if i not in c_ord:
            c_ord.append(i)
    c_ord.reverse()
    
    im_ord = im[c_ord,:,:]
    
    
    # load mask
    mk=imread(mask)
    
    
    # Code from Heath:
    # convert cell mask raster data to polygons for transformation
    cell_shapes = []
    cell_indices = []
    for cell_idx in np.unique(mk)[1:100]:
        cell_mask_thresh = mk.copy()
        cell_mask_thresh[cell_mask_thresh < cell_idx] = 0
        cell_mask_thresh[cell_mask_thresh > cell_idx] = 0
        cell_mask_thresh[cell_mask_thresh == cell_idx] = 255
    
        cell_poly, _ = cv2.findContours(
            cell_mask_thresh.astype(np.uint8), cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE
        )
        if len(cell_poly) == 1:
            cell_shapes.append(np.squeeze(cell_poly).astype(np.double))
            cell_indices.append(cell_idx)
    
    
    # rotate and transform image so that image and mask match
    im_ord = np.rot90(im_ord,axes=[1,2])
    im_ord = np.flip(im_ord, 1)
    
    # add to napari
    # instantiate napari viewer
    viewer = napari.Viewer()
    
    # add image
    viewer.add_image(
        im_ord, channel_axis=0, name=np.array(channels)[c_ord], visible=False, scale=[1, 1]
    )
    # add mask as shape layer
    viewer.add_shapes(cell_shapes,  shape_type="polygon", face_color='#ffffff00', edge_color="#ffffff")
    
    napari.run()
    
if __name__ == '__main__':
    visualize()
