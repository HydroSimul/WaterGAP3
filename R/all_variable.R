#' variable define
#' @name all_vari
#' @description In this topic all of the variable in the EDCHM will be defined.
#'
#'
#' @param n_time (int) number of time steps
#' @param n_spat (int) number of spatial units
#' @param time_dayOfYear_ <1, 366> (d) the number of the day in the year between 1 (1 January) and 365 or 366 (31 December)
#' @param atmos_precipitation_mm (mm/m2/TS) precipitaion volum
#' @param atmos_rain_mm (mm/m2/TS) precipitation in rain form
#' @param atmos_snow_mm (mm/m2/TS) precipitation in snow form
#' @param atmos_temperature_Cel (Cel) the average air temperature in the time phase
#' @param atmos_temperatureMax_Cel (Cel) the maximal air temperature in the time phase
#' @param atmos_temperatureMin_Cel (Cel) the minimal air temperature in the time phase
#' @param atmos_solarRadiat_MJ (MJ/m2/TS) the solar radiation that actually reaches the earths surface
#' @param atmos_netRadiat_MJ (MJ/m2/TS) the balance between the energy absorbed, reflected and emitted by the earths surface or the difference between the incoming net shortwave (Rns) and the net outgoing longwave (Rnl) radiation
#' @param atmos_relativeHumidity_1 (0, 1) relative humidity
#' @param atmos_vaporPress_hPa (hPa) actual vapour pressure, can be calculated by [atmos_VaporPress()]
#' @param atmos_saturatVaporPress_hPa (hPa) saturation vapour pressure at `atmos_temperature_Cel`, can be calculated by [atmos_SaturatVaporPress()]
#' @param atmos_windSpeed_m_s (m/s) measured wind speed at z m above ground surface
#' @param atmos_windMeasureHeight_m (m) height of measurement above ground surface
#' @param atmos_windSpeed2m_m_s (m/s) wind speed at 2 m above ground surface
#' @param atmos_potentialEvatrans_mm (mm/m2/TS) **potential / reference** evapotranspiration
#' @param land_albedo_1 <0, 1> albedo of the region
#' @param land_latitude_Degree (degree) average latitude
#' @param land_elevation_m (m) average elevation
#' @param land_impermeableFrac_1 <0, 1> the maximum impermeable fraction when th soil is fully saturated
#' @param land_water_mm (mm/m2) **pounded water** volume in `landLy` and there is no limit, different than `land_interceptWater_mm`
#' @param land_interceptWater_mm (mm/m2) initial water volume that can be **intercepted**
#' @param land_interceptCapacity_mm (mm/m2) average intercept Capacity (maximal storage capacity)
#' @param land_actualEvatrans_mm (mm/m2/TS) **actual** evapotranspiration from `landLy`
#' @param land_infiltrat_mm (mm/m2/TS) infiltration from `landLy` to `soilLy`
#' @param snow_ice_mm (mm/m2) water equivalent of **ice** in snowpack
#' @param soil_fieldCapacityPerc_1 <0, 1> the relative ratio that the water content can drainage by gravity
#' @param soil_water_mm (mm/m2) water volume in `soilLy`
#' @param soil_capacity_mm (mm/m2) average soil Capacity (maximal storage capacity)
#' @param soil_interflow_mm (mm/m2/TS) subsurface flow directly to the river
#' @param soil_actualEvatrans_mm (mm/m2/TS) **actual** evapotranspiration from `soilLy`
#' @param soil_potentialPercola_mm <0.01, 7> (mm/m2/TS) **potential** percolation
#' @param soil_potentialInteflow_mm <0.01, 7> (mm/m2/TS) **potential** interflow
#' @param soil_potentialCapirise_mm <0.01, 7> (mm/m2/TS) **potential** capillary rise
#' @param ground_water_mm (mm/m2/TS) water volume in `groundLy`
#' @param ground_capacity_mm (mm/m2) water storage capacity in `groundLy`
#' @param ground_lateral_mm (mm/m2/TS) lateral flow, exchange with outside region. It can be **NEGATIV**
#' @param ground_capillarise_mm (mm/m2/TS) capillary rise from `groundLy` to `soilLy`
#' @param ground_potentialLateral_mm <-7, 7> (mm/m2/TS) **potential** lateral flow
#' @param ground_potentialBaseflow_mm <0.01, 7> (mm/m2/TS) **potential** baseflow
#' @param confluen_cellInflow_m3 (m3 / TS) input water volum in every routeline
#' @param confluen_inputWater_mm,,land_runoff_mm,ground_baseflow_mm (mm/m2) input water volum in every routeline
#' @param confluen_iuh_1,confluen_iuhLand_1,confluen_iuhSoil_1,confluen_iuhGround_1 (vector of num, sume() = 1) the ratio in every timestep, can be calculated by [confluenIUH_GR4J1()], [confluenIUH_GR4J2()]
#' @param confluen_responseTime_TS (TS) response or concentration time in every routeline
#' @param water_mm (mm/m2/TS) water volume in `soilLy`, `groundLy` or intercept of `landLy`, same as `soil_water_mm`, `ground_water_mm` or `land_interceptWater_mm`
#' @param capacity_mm (mm/m2) water storage capacity in `soilLy`, `groundLy` or intercept of `landLy`, same as `soil_capacity_mm`, `ground_capacity_mm` or`land_interceptCapacity_mm`
#' @param lake_cellNumber_int (int) lake cell number
#' @param lake_water_m3 (m3) lake water volume
#' @param lake_evatrans_mm (mm) evapotranspiration from the lake
#' @param lake_capacity_m3 (m3) capacity of the lake
#' @param lake_area_km2 (km2) area of the lake
#' @param lake_inflow_m3 (m3/TS) inflow of the lake
#' @param river_water_m3 (m3) river water volume
#' @param river_inflow_m3 (m3/TS) river inflow volume
#' @param river_velocity_km (km/TS) river velocity
#' @param river_length_km (km) river length
#' @param riverlake_cellNumber_int (int) river-lake cell number
#' @param riverlake_area_km2 (km2) area of the river-lake system
#' @param riverlake_capacity_m3 (m3) total capacity of the river-lake system
#' @param riverlake_inflow_m3 (m3) inflow volume to the river-lake system
#' @param riverlake_water_m3 (m3) water volume in the river-lake system
#' @param reservoir_cellNumber_int (int) reservoir cell number
#' @param reservoir_water_m3 (m3) volume of water in the reservoir
#' @param reservoir_capacity_m3 (m3) reservoir capacity
#' @param reservoir_evatrans_mm (mm) evapotranspiration from the reservoir
#' @param reservoir_inflow_m3 (m3/TS) inflow to the reservoir
#' @param reservoir_demand_m3 (m3/TS) water demand from the reservoir
#' @param reservoir_potentialVolume_m3 (m3) potential capacity of the reservoir
#' @param reservoir_area_km2 (km2) surface area of the reservoir
#' @param reservoir_yearInflow_m3 (m3/TS) annual mean inflow to the reservoir
#' @param reservoir_yearDemand_m3 (m3/TS) annual mean water demand from the reservoir
#' @param reservoir_yearRelase_m3 (m3/TS) reference to annual mean water release from the reservoir
#' @param reservoir_isOperateStart_01 (01) indicates if reservoir operation has started (0 or 1)
#' @param reservoir_isIrrigate_01 (01) indicates if the reservoir is used for irrigation (0 or 1)
#' @param basin_landArea_km2 (km2) land area of the basin
#' @param basin_cellNumberStep_int (int) step number of basin cells
#' @param basin_inflowCellNumberStep_int (int) step number of inflow cells in the basin
NULL



