
import os
import numpy as np

from city_metrix.metrix_model import GeoExtent 

import importlib, utils.download as download
importlib.reload(download)

from utils.download import *
from utils.grid import create_grid_for_city

import os, glob
import geopandas as gpd
import pandas as pd




def get_data(city, boundary_path, year, output_base="."):
    # city = "ZAF-Cape_Town"
    # boundary_path = "https://wri-cities-indicators.s3.us-east-1.amazonaws.com/data/published/layers/UrbanExtents/geojson/ZAF-Cape_Town__business_district__UrbanExtents__StartYear_2020_EndYear_2020.geojson"
    data_path = os.path.join(output_base, "data")
    # copy_to_s3 = True

    

    # Get city_polygon -------------------------------------------------
    # city_polygon = get_city_polygon(city, data_path=data_path, copy_to_s3=True)
    city_polygon = gpd.read_file(boundary_path).to_crs(crs='EPSG:4326')

    # Get UTM -------------------------------------------------
    from utils.utm import get_utm
    utm_info = get_utm(city_polygon)
    print(f"UTM EPSG: {utm_info['epsg']}, Earth Engine: {utm_info['ee']}")
    #crs = utm_info['ee']
    crs = 'EPSG:4326'
    
    # Download vector data
    minx, miny, maxx, maxy = city_polygon.total_bounds  
    bbox = GeoExtent(bbox=[minx, miny, maxx, maxy])

    # Create the grid for the city -------------------------------------------------
    city_grid = create_grid_for_city(city, city_polygon, data_path=data_path, copy_to_s3=True)

    # Start a dask client
    from dask.distributed import Client
    from dask import delayed
    import dask
    client = Client()

    # Print Dask client information
    print(f"Dask: {client}")
    print(f"Dask dashboard: {client.dashboard_link}")

    try:
        # Create all tasks for all cells and all data types to run in parallel
        print("Creating all data fetching tasks for all grid cells...")
        lst_alb_tasks = []
        fr_tasks = []

        for idx, cell in city_grid.iterrows():
            grid_cell_id = cell['ID']
            geometry = cell['geometry']
            # Get the bounding box of the cell
            bbox = GeoExtent(bbox=geometry.bounds)
            
            geographic_box = bbox.as_geographic_bbox()
            geographic_bbox_str = ','.join(map(str, geographic_box.bounds))
            print(f"Preparing tasks for grid cell {grid_cell_id} with bbox: {geographic_bbox_str}")
            
            # Create all data fetching tasks for this cell
            lst_alb_tasks.append(delayed(get_lst)(
                city,
                geographic_box,
                grid_cell_id=grid_cell_id,
                data_path=data_path,
                center_date="2022-01-22",
                copy_to_s3=False,
            ))
            lst_alb_tasks.append(delayed(get_albedo)(
                city,
                geographic_box,
                grid_cell_id=grid_cell_id,
                data_path=data_path,
                year=year,
                copy_to_s3=False,
            ))
        
            # Fractional vegetation tasks are delayed but stored separately
            fr_tasks.append(delayed(get_fr)(
                city,
                geographic_box,
                year=year,
                grid_cell_id=grid_cell_id,
                data_path=data_path,
                copy_to_s3=False,
            ))



            # Add all tasks to the master list
            # all_tasks.extend([lst_task, fr_task, alb_task])

            # Keep track of what each task does for debugging
            # task_descriptions.extend([
            #     f"lst_cell_{grid_cell_id}",
            #     f"fr_cell_{grid_cell_id}",
            #     f"alb_cell_{grid_cell_id}"
            # ])
        
        # Execute ALL tasks in parallel across all cells and all data types
        # print(f"Executing ALL {len(all_tasks)} tasks in parallel...")
        # print(f"This includes {len(city_grid)} cells × 6 data types = {len(all_tasks)} total tasks")
        # print("=" * 60)

        visualize = False
        if visualize:
            print("Visualizing task graph...")
            try:
                dask.visualize(*all_tasks, filename='task_graph.png', optimize_graph=False)
                print("Task graph saved as 'task_graph.png' - check it out!")
            except Exception as e:
                print(f"Visualization failed: {e}")
        else:
        
            # Now execute the actual computation
            print("Starting actual computation...")
            results = dask.compute(*lst_alb_tasks)
            print("Computing fractional vegetation sequentially…")
            for fr in fr_tasks:
                fr.compute()

            print("ALL data fetching completed!")
            print(f"Successfully executed {len(results)} tasks in parallel")
            print("=" * 60)

    finally:
        # Close the Dask client
        client.close()
        
    

# if __name__ == '__main__':
#     main()

