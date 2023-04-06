
library(tidyverse)
library(sf)
library(fable)
library(fabletools)
library(feasts)
library(data.table)
library(terra)
library(lubridate)
library(tsibble)
library(ggplot2)
library(scales)
library(zoo)
library(plotly)
library(shiny)
library(leaflet)


model_vegetation <- function(refl_l8, gpm){
  
  refl_l8$date <- as.Date(refl_l8$date, tryFormats = "%d/%m/%Y")
  gpm$date <- as.Date(gpm$date, tryFormats = "%d/%m/%Y")
  
  
  # Pre-process
  cloud_landsat <- function(x){
    bits <- as.numeric(intToBits(x))[1:16]
    if(bits[7] == 1 && bits[5] == 0 && bits[3] == 0){
      return(0)
    }
    else(return(1))
  }
  
  ######## Convert reflectances and gpm values
  # CONVERT LANDSAT 8
  # Scale DN values to reflectance values
  scale_fac <- 0.0000275
  offset <- -0.2
  refl_l8[,c('SR_B4', 'SR_B5')] <- (refl_l8[,c('SR_B4', 'SR_B5')] * scale_fac) + offset
  
  # Remove negative reflectance
  refl_l8 <- refl_l8[refl_l8$SR_B4 >= 0,]
  refl_l8 <- refl_l8[refl_l8$SR_B5 >= 0,]
  
  # Remove clouds
  refl_l8$clouds <- 0
  refl_l8$clouds <- lapply(refl_l8$QA_PIXEL, cloud_landsat) # Apply function to all pixels 
  refl_l8 <- refl_l8[refl_l8$clouds==0,] # Filter out pixels with clouds
  
  # Compute NDVI
  refl_l8$ndvi <- (refl_l8$SR_B5 - refl_l8$SR_B4) / (refl_l8$SR_B5 + refl_l8$SR_B4)
  
  # CONVERT GPM
  gpm$precipitationCal <- 0.5 * gpm$precipitationCal
  gpm <- gpm %>%   
    mutate(yearmonth = yearmonth(date))
  
  # Merge GPM and Landsat data
  ts_month <- refl_l8 %>%
    group_by(yearmonth(date)) %>% 
    summarise(ndvi = mean(ndvi, na.rm=T)) %>% 
    rename(yearmonth = `yearmonth(date)`) %>% 
    right_join(gpm[,c('precipitationCal', 'yearmonth')], by='yearmonth') %>% 
    rename(prcp = precipitationCal) %>% 
    arrange(yearmonth)
  
  
  ############## Creates tsibble with lagged precipitation and interpolate the missing ndvi values
  ts_month_int <- ts_month %>% 
    # Add lagged precipitation variables
    mutate(prcp_lag1=data.table::shift(prcp, n=1, type='lag'),
           prcp_lag2=data.table::shift(prcp, n=2, type='lag')) %>% 
    # Check if plot ts is all NAs and exclude
    mutate(entire_na = length(which(!is.na(ndvi))) == 0) %>% 
    filter(!entire_na) %>% 
    # Exclude NAs when at beginning or end of time series
    slice(min(which(!is.na(ndvi))):max(which(!is.na(ndvi)))) %>% 
    # Remove the entire.na column
    select(!entire_na) %>% 
    # Interpolate NDVI values
    mutate(ndvi_int=na.approx(ndvi)) %>% 
    # Convert to tsibble
    as_tsibble(index=yearmonth)
  
  
  
  ############ Train model
  ts_ref <- ts_month_int %>% filter(year(yearmonth) < 2017)
  ts_pred <- ts_month_int %>% filter(year(yearmonth) >= 2017)
  
  armax <- ts_ref %>% 
    model(ARIMA(ndvi_int ~ prcp + prcp_lag1 + prcp_lag2 + pdq(d=0) + PDQ(D=0), stepwise = T, ic='aicc'))
  
  # Forecast ndvi based on models and precipitation
  fc <- fabletools::forecast(armax, new_data = ts_pred[,c('yearmonth', 'prcp', 'prcp_lag1', 'prcp_lag2')])
  
  # Extract 95% confidence levels
  ci <- fc$ndvi_int %>% 
    hilo(80) 
  fc <- fc %>% 
    mutate(upper = ci$upper,
           lower = ci$lower)
  
  
  ################# Add forecast and actual ndvi
  fc_dt <- data.table(fc[,c('.mean', 'yearmonth','upper', 'lower')])
  names(fc_dt)[names(fc_dt) == '.mean'] <- 'pred_test'
  # Merge actual ndvi with the predicted ndvi
  fc_ts <- merge(ts_month_int, fc_dt, by=c('yearmonth'), all.x=T)
  # Merge actual ndvi with predicted reference ndvi
  fc_ts <- merge(fc_ts, residuals(armax)[c('yearmonth', '.resid')], 
                 by=c('yearmonth'), all.x=T)
  # Compute training prediction
  fc_ts$pred_train <- fc_ts$ndvi_int + fc_ts$.resid
  # Combine test and training prediction into one and convert to tsibble
  fc_ts <- fc_ts %>% 
    tsibble(index=yearmonth) %>% 
    mutate(ndvi_pred = if_else(is.na(pred_test), pred_train, pred_test))
  
  
  ############### Compute MAE
  mae <- fc_ts %>% 
    mutate(test_res = if_else(ndvi_int > upper, ndvi_int - upper, 0),
           training = if_else(!is.na(pred_train), 1, 0)) %>% 
    tibble() %>% 
    summarise(n = n(),
              n_training = sum(training),
              mae = sum(test_res, na.rm=T) / (n - n_training))
  
  
  ########## Plot
  fc_ci <- fc %>% 
    mutate(ci_95 = hilo(ndvi_int,95),
           upper_95=ci_95$upper,
           lower_95=ci_95$lower,
           ci_80 = hilo(ndvi_int, 80),
           upper_80=ci_80$upper,
           lower_80=ci_80$lower,
           ci_50 = hilo(ndvi_int, 50),
           upper_50=ci_50$upper,
           lower_50=ci_50$lower) 
  
  restoration_plot <- plot_ly(fc_ts, x = ~as.Date(yearmonth), y = ~ndvi_int, name = 'NDVI', type = 'scatter', mode = 'lines',line = 
                                list(color='rgb(0,100,80)'), alpha=1) %>% 
    add_lines(x = ~as.Date(yearmonth), y = ~ndvi_pred, name = 'NDVI prediction', type = 'scatter', mode = 'lines',line
              =list(color='rgb(0,0,60)'),alpha=0.5)%>%
    # 80% CI
    add_ribbons(fc_ts, x = ~as.Date(yearmonth), ymin = ~lower, ymax = ~upper, line = list(color = 'rgba(7, 164, 181, 0)'), 
                fillcolor = 'rgba(7, 164, 181, 0.2)', name = '80% confidence', alpha=0.1) %>% 
    # 95% CI
    add_ribbons(data = fc_ci, x = ~as.Date(yearmonth), ymin = ~lower_95, ymax = ~upper_95, line = list(color = 'rgba(7, 164, 181,0)'), 
                fillcolor = 'rgba(7,164, 181, 0.2)', name = '95% confidence',alpha=0.2) %>% 
    add_annotations(x=as.Date('2022-01-01'), y=0.7, text=paste0('MAEp = ', round(mae,3)), showarrow=F) %>% 
    layout(yaxis = list(range=c(0,1), title = 'NDVI'),
           xaxis = list(title = ''))
  
  # Return plot
  return(restoration_plot)
}

