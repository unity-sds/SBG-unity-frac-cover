#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Winston Olson-Duvall
"""
import datetime as dt
import glob
import json
import os
import subprocess
import shutil
import sys
import hytools as ht
import pandas as pd
import numpy as np
import pystac
from pathlib import Path

try:
    from osgeo import gdal
except:
    import gdal
from PIL import Image

# stage_in packages
from unity_sds_client.resources.collection import Collection

# stage_out packages
from datetime import datetime, timezone
from unity_sds_client.resources.dataset import Dataset
from unity_sds_client.resources.data_file import DataFile
import spectral.io.envi as envi


def find_nearest_day(day, subset_days):
    distances = []
    for special_day in subset_days:
        diff = (special_day - day) % 365
        distances.append(min(diff, 365 - diff))
    nearest_day = subset_days[distances.index(min(distances))]
    return np.argwhere(nearest_day==subset_days)[0][0]

def get_frcov_basename(corfl_basename, crid):
    # Replace product type
    tmp_basename = corfl_basename.replace("L2A_CORFL", "L2B_FRCOV")
    # Split, remove old CRID, and add new one
    tokens = tmp_basename.split("_")[:-1] + [str(crid)]
    return "_".join(tokens)

'''
def generate_stac_metadata(basename, description, in_meta):

    out_meta = {}
    out_meta['id'] = basename
    out_meta['start_datetime'] = dt.datetime.strptime(in_meta['start_datetime'], "%Y-%m-%dT%H:%M:%SZ")
    out_meta['end_datetime'] = dt.datetime.strptime(in_meta['end_datetime'], "%Y-%m-%dT%H:%M:%SZ")
    out_meta['geometry'] = in_meta['geometry']
    base_tokens = basename.split('_')
    out_meta['collection'] = f"SISTER_{base_tokens[1]}_{base_tokens[2]}_{base_tokens[3]}_{base_tokens[5]}"
    out_meta['properties'] = {
        'sensor': in_meta['sensor'],
        'description': description,
        'product': base_tokens[3],
        'processing_level': base_tokens[2],
        "cover_percentile_counts": in_meta["cover_percentile_counts"]
    }
    return out_meta
'''

def create_item(metadata, assets):
    item = pystac.Item(
        id=metadata['id'],
        datetime=metadata['start_datetime'],
        start_datetime=metadata['start_datetime'],
        end_datetime=metadata['end_datetime'],
        geometry=metadata['geometry'],
        collection=metadata['collection'],
        bbox=None,
        properties=metadata['properties']
    )
    # Add assets
    for key, href in assets.items():
        item.add_asset(key=key, asset=pystac.Asset(href=href))
    return item



"""
This is the main function

"""

# The defaults used here generally relflect a local or jupyter environment; they are replaced with "runtime" values when run in the system.
input_stac_collection_file = sys.argv[1] #'/unity/ads/input_collections/SBG-L2-FRAC_COVER/catalog.json' # type: stage-in
output_stac_catalog_dir    = sys.argv[2] #'/unity/ads/outputs/SBG-L2-FRAC_COVER/'                    # type: stage-out
n_cores = sys.argv[3] #"1"
refl_scale = sys.argv[4] #"1.0"
normalization = sys.argv[5] #"none"
crid = sys.argv[6] #"001"
experimental_flag = sys.argv[7].lower() == 'true' #sys.argv[7] #"True"
output_collection_name = sys.argv[8] #'SBG-L2-Fractional-Cover'

#temp work dir
#optional variables
temp_work_dir = "/tmp/frcover"
if not os.path.exists(temp_work_dir):
    os.makedirs(temp_work_dir+"/work", exist_ok=True)

# # Import Files from STAC Item Collection
# 
# Load filenames from the stage_in STAC item collection file
inp_collection = Collection.from_stac(input_stac_collection_file)
data_filenames = inp_collection.data_locations(["data"])

print(f"Data Files (JSON): {data_filenames}")

# Define paths and variables
from pathlib import Path
sister_frcov_dir = os.path.abspath(os.path.dirname(__file__))
print(sister_frcov_dir)
# Make work dir
# if not os.path.exists(f"{sister_frcov_dir}/work"):
#     print("Making work directory")
#     os.mkdir(f"{sister_frcov_dir}/work")
        
if experimental_flag:
    disclaimer = "(DISCLAIMER: THIS DATA IS EXPERIMENTAL AND NOT INTENDED FOR SCIENTIFIC USE) "
else:
    disclaimer = ""


specun_dir = os.path.join(os.path.dirname(sister_frcov_dir), "SpectralUnmixing")
print(f'Specun_dir: {specun_dir}')
corfl_file_path = os.path.dirname(data_filenames[0])
corfl_basename = Path(data_filenames[0]).stem
frcov_basename = get_frcov_basename(corfl_basename, crid)
corfl_img_path = f"{temp_work_dir}/work/{corfl_basename}"
corfl_hdr_path = f"{temp_work_dir}/work/{corfl_basename}.hdr"
frcov_img_path = f"{temp_work_dir}/work/{frcov_basename}"

# Copy the input files into the work directory (don't use .bin)
shutil.copyfile(f"{corfl_file_path}/{corfl_basename}.bin", corfl_img_path)
shutil.copyfile(f"{corfl_file_path}/{corfl_basename}.hdr", corfl_hdr_path)


#Load reflectance im
rfl = ht.HyTools()
rfl.read_file(corfl_img_path)

line_data = (rfl.get_band(0) == rfl.no_data).sum(axis=1)
start_line = 1+np.argwhere(line_data != rfl.columns)[0][0]
end_line = rfl.lines - np.argwhere(np.flip(line_data) != rfl.columns)[0][0] -1

endmember_lib_path = f"{sister_frcov_dir}/data/veg_soil_water_snow_endmembers.csv"
print(f'endmember: {endmember_lib_path} ')
endmembers = pd.read_csv(endmember_lib_path)

no_snow = endmembers[endmembers['class'] != 'snow']
no_snow.to_csv(f'{temp_work_dir}/work/endmembers_no_snow.csv',index = False)

no_water = endmembers[endmembers['class'] != 'water']
no_water.to_csv(f'{temp_work_dir}/work/endmembers_no_water.csv',index = False)

# Build command and run unmix.jl
unmix_exe = f"{specun_dir}/unmix.jl"
log_path = f"{output_stac_catalog_dir}/{frcov_basename}.log"

snow_clim_file = f'{sister_frcov_dir}/data/LIN10A1_snow_climatology_13day.tif'
snow_clim = gdal.Open(snow_clim_file)
snow = snow_clim.GetRasterBand(1)
snow_days = np.arange(1,365,13)

ulx,pixel_size,_,uly,_,_ = snow_clim.GetGeoTransform()

run_config = {}
stac_json_path = os.path.join(corfl_file_path, f"{corfl_basename}.json")
with open(stac_json_path, "r") as f:
    stac_item = json.load(f)
run_config["metadata"] = stac_item["properties"]
run_config["metadata"]["geometry"] = stac_item["geometry"]

bbox_min_x,bbox_min_y = np.min(run_config['metadata']['geometry']['coordinates'][:4],axis=0)
bbox_max_x,bbox_max_y = np.max(run_config['metadata']['geometry']['coordinates'][:4],axis=0)

x_offset = int((bbox_min_x-ulx)/pixel_size)
width = int((bbox_max_x-ulx)/pixel_size) -x_offset

y_offset = int((uly-bbox_max_y)/pixel_size)
height = int((uly-bbox_min_y)/pixel_size) - y_offset

subset = snow.ReadAsArray(x_offset, y_offset,
                          width, height)

datetime = dt.datetime.strptime(run_config['metadata']['start_datetime'], '%Y-%m-%dT%H:%M:%SZ')
doy = datetime.timetuple().tm_yday
bit = find_nearest_day(doy, snow_days)
snow_mask =(subset >> bit) &1
snow_present = snow_mask.sum() >0

if snow_present:
    print('Snow expected')
    endmember_files = ['no_water','no_snow']
else:
    endmember_files = ['no_snow']

for endmember_file in endmember_files:
    # Add required args
    cmd = [f'JULIA_PROJECT={specun_dir}', "julia"]

    # Add parallelization if n_cores > 1
    if int(n_cores) > 1:
        cmd += ["-p", n_cores]

    cmd += [
        unmix_exe,
        corfl_img_path,
        f'{temp_work_dir}/work/endmembers_{endmember_file}.csv',
        "class",
        f"{frcov_img_path}_{endmember_file}",
        "--mode=sma",
        f"--log_file={log_path}",
        "--num_endmembers=3",
        f"--refl_nodata={rfl.no_data}",
        f"--start_line={start_line}",
        f"--end_line={end_line}",
        "--normalization=brightness",
        f"--refl_scale={refl_scale}"]

    print("Running unmix.jl command: " + " ".join(cmd))
    import subprocess, os
    my_env = os.environ.copy()
    subprocess.run(" ".join(cmd), shell=True, env=my_env)

## Unity (BLEE): Modified due to lack of files:
no_snow_frcov_file = f'{frcov_img_path}_no_snow_fractional_cover'
if os.path.isfile(no_snow_frcov_file):
    no_snow_gdal = gdal.Open(no_snow_frcov_file)
    no_snow_frcov  = no_snow_gdal.ReadAsArray()
    no_snow_frcov[no_snow_frcov==rfl.no_data] = np.nan
else:
    print(f'WARNING: Following file could not be found: {no_snow_frcov_file}.')

filter_frcov = np.zeros((rfl.lines,rfl.columns,4))

no_water_frcov_file = f'{frcov_img_path}_no_water_fractional_cover'
if snow_present:
    if os.path.isfile(no_water_frcov_file):
        no_water_gdal = gdal.Open(no_water_frcov_file)
        no_water_frcov  = no_water_gdal.ReadAsArray()
        no_water_frcov[no_water_frcov==rfl.no_data] = np.nan

        water = (rfl.ndi(550,850) > 0) & (rfl.get_wave(900) < .15)

        filter_frcov[water,:3] =no_snow_frcov[:3,water].T
        filter_frcov[~water,0] =no_water_frcov[0,~water].T
        filter_frcov[~water,1] =no_water_frcov[1,~water].T
        filter_frcov[~water,3] =no_water_frcov[2,~water].T
        
    else: 
        print(f'WARNING: Following file could not be found: {no_water_frcov_file}.')

elif os.path.isfile(no_snow_frcov_file):
    filter_frcov[:,:,:3]  = np.moveaxis(no_snow_frcov[:3],0,-1)
##
    
filter_frcov[~rfl.mask['no_data']] = -9999
filter_frcov_file = f'{frcov_img_path}_fractional_cover'
band_names = ['soil','vegetation','water','snow_ice']
header = rfl.get_header()
header['bands'] = 4
header['band names'] = band_names
header['wavelength'] = []
header['fwhm'] = []
header['bbl'] = []

writer = ht.io.envi.WriteENVI(filter_frcov_file, header)
for band_num in range(4):
    writer.write_band(filter_frcov[:,:,band_num], band_num)

cover_counts = {}

# Load bands and calculate cover percentile counts
for band_num in range(0,4):
    frac_dict = {}
    if ("DESIS" in corfl_img_path) and (band_num ==3):
        band_arr = np.zeros((rfl.lines,rfl.columns))
    else:
        band_arr = filter_frcov[:,:,band_num]
    for percent in np.linspace(0,1,11):
        frac_dict[round(percent,1)] = float((band_arr >= percent).sum())
    cover_counts[band_names[band_num]] = frac_dict

run_config["metadata"]["cover_percentile_counts"] = cover_counts

# Generate quicklook
frcov_ql_path = f"{output_stac_catalog_dir}/{frcov_basename}.png"
print(f"Generating quicklook to {frcov_ql_path}")

rgb=  np.array(filter_frcov[:,:,:3])
rgb[rgb == rfl.no_data] = np.nan
rgb = (rgb*255).astype(np.uint8)
im = Image.fromarray(rgb)
im.save(frcov_ql_path)

#Convert to COG
temp_file =  f'{temp_work_dir}/work/temp_frcover.tif'
out_file =  f"{output_stac_catalog_dir}/{frcov_basename}.tif"

print(f"Creating COG {out_file}")

units = ['PERCENT',
        'PERCENT',
        'PERCENT',
        'PERCENT',]

descriptions=  ['SOIL PERCENT COVER',
                  'VEGETATION PERCENT COVER',
                  'WATER PERCENTCOVER',
                  'SNOW/ICE PERCENT COVER']

# Set the output raster transform and projection properties
driver = gdal.GetDriverByName("GTIFF")
tiff = driver.Create(temp_file,
                      no_snow_gdal.RasterXSize,
                      no_snow_gdal.RasterYSize,
                      4,
                      gdal.GDT_Float32)

tiff.SetGeoTransform(no_snow_gdal.GetGeoTransform())
tiff.SetProjection(no_snow_gdal.GetProjection())
fc_description = f"{disclaimer}FRACTIONAL COVER"
tiff.SetMetadataItem("DESCRIPTION", fc_description)

# Write bands to file
for i,band_name in enumerate(band_names):
    out_band = tiff.GetRasterBand(i+1)
    out_band.WriteArray(filter_frcov[:,:,i])
    out_band.SetDescription(band_name)
    out_band.SetNoDataValue(rfl.no_data)
    out_band.SetMetadataItem("UNITS",units[i])
    out_band.SetMetadataItem("DESCRIPTION",descriptions[i])
del tiff, driver

os.system(f"gdaladdo -minsize 900 {temp_file}")
os.system(f"gdal_translate {temp_file} {out_file} -co COMPRESS=LZW -co TILED=YES -co COPY_SRC_OVERVIEWS=YES")

# If experimental, prefix filenames with "EXPERIMENTAL-"
if experimental_flag:
    for file in glob.glob(f"{output_stac_catalog_dir}/SISTER*"):
        shutil.move(file, f"{output_stac_catalog_dir}/EXPERIMENTAL-{os.path.basename(file)}")

# Update the path variables if now experimental
frcov_file = glob.glob(f"{output_stac_catalog_dir}/*{crid}.tif")[0]
#out_runconfig = glob.glob(f"{output_stac_catalog_dir}/*{crid}.runconfig.json")[0]
log_path = glob.glob(f"{output_stac_catalog_dir}/*{crid}.log")[0]
frcov_basename = os.path.basename(frcov_file)[:-4]


'''
## Removed for Unity. Is this needed? 
# Generate STAC
catalog = pystac.Catalog(id=frcov_basename,
                         description=f'{disclaimer}This catalog contains the output data products of the SISTER '
                                     f'fractional cover PGE, including a fractional cover cloud-optimized GeoTIFF. '
                                     f'Execution artifacts including the runconfig file and execution '
                                     f'log file are included with the fractional cover data.')

# Add items for data products
tif_files = glob.glob(f"{output_stac_catalog_dir}/*SISTER*.tif")
tif_files.sort()
for tif_file in tif_files:
    metadata = generate_stac_metadata(frcov_basename, fc_description, run_config["metadata"])
    assets = {
        "cog": f"./{os.path.basename(tif_file)}",
    }
    # If it's the fractional cover product, then add png, runconfig, and log
    if os.path.basename(tif_file) == f"{frcov_basename}.tif":
        png_file = tif_file.replace(".tif", ".png")
        assets["browse"] = f"./{os.path.basename(png_file)}"
        #assets["runconfig"] = f"./{os.path.basename(out_runconfig)}"
        if os.path.exists(log_path):
            assets["log"] = f"./{os.path.basename(log_path)}"
    item = create_item(metadata, assets)
    catalog.add_item(item)

# set catalog hrefs
catalog.normalize_hrefs(f"./{output_stac_catalog_dir}/{frcov_basename}")

# save the catalog
catalog.describe()
catalog.save(catalog_type=pystac.CatalogType.SELF_CONTAINED)
print("Catalog HREF: ", catalog.get_self_href())

# Move the assets from the output directory to the stac item directories and create empty .met.json files
for item in catalog.get_items():
    for asset in item.assets.values():
        fname = os.path.basename(asset.href)
        shutil.move(f"{output_stac_catalog_dir}/{fname}", f"{output_stac_catalog_dir}/{frcov_basename}/{item.id}/{fname}")
    with open(f"{output_stac_catalog_dir}/{frcov_basename}/{item.id}/{item.id}.met.json", mode="w"):
        pass
'''   

## Create stage-out item catalog

#Uncertain about this part. What is supposed to get sent out to STAC? 
from datetime import datetime, timezone

# Create a collection
out_collection = Collection(output_collection_name)

# Add output file(s) to the dataset
file = glob.glob(f"{output_stac_catalog_dir}/*{crid}*")

if file:   
    start_time = dt.datetime.strptime(run_config["metadata"]['start_datetime'], "%Y-%m-%dt%H:%M:%Sz")
    end_time = dt.datetime.strptime(run_config["metadata"]['end_datetime'], "%Y-%m-%dt%H:%M:%Sz")
    name = os.path.splitext(os.path.basename(file[0]))[0]
    description = f"{disclaimer}This catalog contains the output data products of the SISTER fractional cover PGE, including a fractional cover cloud-optimized GeoTIFF. Execution artifacts including the runconfig file and execution log file are included with the fractional cover data."

    dataset = Dataset(
        name=name,
        collection_id=out_collection.collection_id, 
        start_time=start_time.strftime("%Y-%m-%dT%H:%M:%SZ"),
        end_time=end_time.strftime("%Y-%m-%dT%H:%M:%SZ"),
        creation_time=datetime.utcnow().replace(tzinfo=timezone.utc).isoformat(),
        )

    for file in glob.glob(f"{output_stac_catalog_dir}/*{crid}*"):  

        if file.endswith(".bin"):
            dataset.add_data_file(DataFile("binary", file, ["data"]))
        elif file.endswith(".png"):
            dataset.add_data_file(DataFile("image/png", file, ["browse"]))
        elif file.endswith(".tif"):
            dataset.add_data_file(DataFile("image/tif", file, ["data"]))
        elif file.endswith(".hdr"):
            dataset.add_data_file(DataFile("header", file, ["data"]))
        else:
            dataset.add_data_file(DataFile(None, file, ["metadata"]))

    dataset.add_data_file(DataFile("text/json", output_stac_catalog_dir + '/' +  name +'.json', ["metadata"]))
    dataset.add_property("Description", description)

    for k in run_config["metadata"]["cover_percentile_counts"].keys():
        dataset.add_property(k, run_config["metadata"]["cover_percentile_counts"][k])

# Add the dataset to the collection
out_collection._datasets.append(dataset)

Collection.to_stac(out_collection, output_stac_catalog_dir)

