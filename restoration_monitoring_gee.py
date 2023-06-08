# %% 
# Load packages
  
import ee
import pandas as pd
from datetime import datetime, date, time, timedelta
from dateutil.relativedelta import *
import geopandas as gpd
import numpy as np
import datatable as dt
from datatable import f
import time
from pprint import pprint
import os, re
from gee_subset import gee_subset
import shapely.geometry
import re
import geemap
from shapely.geometry import Polygon

# Trigger the authentication flow.
# ee.Authenticate()
# Initialize the library.
ee.Initialize()


# %% 
# Put everything in a function
def extract_gee(lon, lat):

  # Geometry
  # Store in datatable
  plots_dt = dt.Frame({"lon": [lon], "lat": [lat]})

  point = [plots_dt[0, 'lon'], plots_dt[0, 'lat']]
  centroid_multi = ee.Geometry.MultiPoint(point)

  start_date = datetime(2013, 1, 1)
  end_date = datetime.today()

  # Select bands and resolution to extract from landsat 8
  bands = ["SR_B4", "SR_B5", "QA_PIXEL"]
  scale = 30
  product='LANDSAT/LC08/C02/T1_L2'

  # EXTRACT LANDSAT 8
  geometry = ee.Geometry.Point([plots_dt[0,'lon'], plots_dt[0, 'lat']])
  col = ee.ImageCollection(product).\
        select(tuple(bands)).\
        filterDate(start_date, end_date).filterBounds(geometry)
      
  region = col.getRegion(geometry, int(scale)).getInfo()
  df = pd.DataFrame.from_records(region[1:len(region)])
  df.columns = region[0]
  df = df[['time', 'SR_B4', 'SR_B5', 'QA_PIXEL']]
  
  df.time = df.time / 1000
  df['time'] = pd.to_datetime(df['time'], unit = 's')
  df.rename(columns = {'time': 'date'}, inplace = True)
  df.sort_values(by = 'date')

  # Transform to dt
  l8_out = dt.Frame(df)


  # EXTRACT GPM
  # Specify number of days in period of interest
  d0 = datetime(start_date.year, start_date.month, start_date.day)
  d1 = datetime(end_date.year, end_date.month, 1)
  delta = d1 - d0
  days = delta.days

  # number of months in period
  def diff_month(d1, d2):
     return (d1.year - d2.year) * 12 + d1.month - d2.month
  months_ts = diff_month(d1, d0)

  # Create list with the dates off all the observations
  months_date = []
  for m in range(months_ts):
    first_month = start_date
    next_month = first_month + relativedelta(months =+ m)
    months_date.append(next_month)

  # Extract all the precipitation data for all sites
  gpm = ee.ImageCollection('NASA/GPM_L3/IMERG_V06').\
         select('precipitationCal').\
         filterDate(start_date, end_date).filterBounds(centroid_multi)

  # Create a function to go over the FeatureCollection and take the monthly sum
  def GPMsum(img_collection):
    mylist = ee.List([])
    for m in range(months_ts):

      ini = start_date + relativedelta(months=+m)
      end = ini + relativedelta(months=+1) + relativedelta(days=-1)

      sum_image = img_collection.filterDate(ini,end).select(0).sum()
      mylist = mylist.add(sum_image.set('system:time_start', ini))
    return ee.ImageCollection.fromImages(mylist)

  # Apply the 'GPMsum' function to create FeatureCollection with monthly sum
  monthlyGPM = ee.ImageCollection(GPMsum(gpm))

  # Extract the GPM values
  region = monthlyGPM.getRegion(centroid_multi, int(scale)).getInfo()
  df = pd.DataFrame.from_records(region[1:len(region)])
  df.columns = region[0]
  df = df[['id', 'precipitationCal']]
  df = dt.Frame(df)
  # Add datetime
  df['date'] = dt.Frame(months_date)

  return l8_out, df

#####################################################################################
#####################################################################################
#####################################################################################


# Extract GEE for polygon

def extract_gee_polygon(cor1, cor2, cor3, cor4, cor5):

  ### Make polygon geometry
  cords = [[cor1,cor2,cor3,cor4,cor5]]
  geometry = ee.Geometry.Polygon(cords)
  
  ### Make point geometry
  # Compute centroid longitude and latitude
  lon = np.mean(np.array([cor1[0], cor3[0]]))
  lat = np.mean(np.array([cor1[1], cor3[1]]))
  
  # Store in datatable
  plots_dt = dt.Frame({"lon": [lon], "lat": [lat]})
  
  # Make ee.Point
  point = [plots_dt[0, 'lon'], plots_dt[0, 'lat']]
  centroid_multi = ee.Geometry.MultiPoint(point)

  start_date = datetime(2013, 1, 1)
  end_date = datetime.today()

  # Select bands and resolution to extract from landsat 8
  bands = ["SR_B4", "SR_B5", "QA_PIXEL"]
  scale = 30
  product='LANDSAT/LC08/C02/T1_L2'

  # 1) EXTRACT LANDSAT 8
  col = ee.ImageCollection(product).\
      select(tuple(bands)).\
      filterDate(start_date, end_date).filterBounds(geometry)

  region = col.getRegion(geometry, int(scale)).getInfo()

  # If no pixels in geometry, take centroid of plot
  if len(region) == 1:
      geometry = ee.Geometry.Point([plots_dt[0,'lon'], plots_dt[0, 'lat']])
      col = ee.ImageCollection(product).\
          select(tuple(bands)).\
          filterDate(start_date, end_date).filterBounds(geometry)
    
  # Create df from ee.ImageCollection  
  region = col.getRegion(geometry, int(scale)).getInfo()
  df = pd.DataFrame.from_records(region[1:len(region)])
  df.columns = region[0]
  df = df[['time', 'SR_B4', 'SR_B5', 'QA_PIXEL']]
  
  df.time = df.time / 1000
  df['time'] = pd.to_datetime(df['time'], unit = 's')
  df.rename(columns = {'time': 'date'}, inplace = True)
  df.sort_values(by = 'date')

  # Transform to dt
  l8_out = dt.Frame(df)


  # 2) EXTRACT GPM
  # Specify number of days in period of interest
  d0 = datetime(start_date.year, start_date.month, start_date.day)
  d1 = datetime(end_date.year, end_date.month, 1)
  delta = d1 - d0
  days = delta.days

  # number of months in period
  def diff_month(d1, d2):
     return (d1.year - d2.year) * 12 + d1.month - d2.month
  months_ts = diff_month(d1, d0)

  # Create list with the dates off all the observations
  months_date = []
  for m in range(months_ts):
    first_month = start_date
    next_month = first_month + relativedelta(months =+ m)
    months_date.append(next_month)

  # Extract all the precipitation data for all sites
  gpm = ee.ImageCollection('NASA/GPM_L3/IMERG_V06').\
         select('precipitationCal').\
         filterDate(start_date, end_date).filterBounds(centroid_multi)

  # Create a function to go over the FeatureCollection and take the monthly sum
  def GPMsum(img_collection):
    mylist = ee.List([])
    for m in range(months_ts):

      ini = start_date + relativedelta(months=+m)
      end = ini + relativedelta(months=+1) + relativedelta(days=-1)

      sum_image = img_collection.filterDate(ini,end).select(0).sum()
      mylist = mylist.add(sum_image.set('system:time_start', ini))
    return ee.ImageCollection.fromImages(mylist)

  # Apply the 'GPMsum' function to create FeatureCollection with monthly sum
  monthlyGPM = ee.ImageCollection(GPMsum(gpm))

  # Extract the GPM values
  region = monthlyGPM.getRegion(centroid_multi, int(scale)).getInfo()
  df = pd.DataFrame.from_records(region[1:len(region)])
  df.columns = region[0]
  df = df[['id', 'precipitationCal']]
  df = dt.Frame(df)
  # Add datetime
  df['date'] = dt.Frame(months_date)

  return l8_out, df