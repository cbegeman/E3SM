import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

run_path = '/lcrc/group/acme/ac.dcomeau/E3SMv3_dev/20240201.v3.SORRM.piControl.chrysalis/run/'
file_prefix = '20240201.v3.SORRM.piControl.chrysalis.'
conserv_prefix = 'mpaso.hist.am.conservationCheck'
ts_prefix = 'mpaso.hist.am.timeSeriesStatsMonthly'
output_dir = '/lcrc/group/acme/ac.cbegeman/ssh_analysis'

decimal_year = []
total_ssh = []
ssh_rate = []
A = []
V = []
accumulated_snow_flux = []
accumulated_rain_flux = []
accumulated_evaporation_flux = []
accumulated_sea_ice_flux = []
accumulated_land_ice_flux = []
accumulated_frazil_flux = []
accumulated_iceberg_flux = []
accumulated_river_runoff_flux = []
accumulated_ice_runoff_flux = []
accumulated_removed_river_runoff_flux = []
accumulated_removed_ice_runoff_flux = []

s_year = 3.154e7 
s_month = s_year / 12.
rho_sw = 1026.
mass_flux_to_SLR = 1.e3 / rho_sw # we still need to divide by area of the ocean below to get units of mm/s
for year in range(30):
   for month in range(12):
       decimal_year.append(year + 1. + month / 12.)
       date_suffix = f'{year+1:04g}-{month+1:02g}-01.nc'
       ds = xr.open_dataset(f'{run_path}/{file_prefix}{conserv_prefix}.{date_suffix}')
       ds_ts = xr.open_dataset(f'{run_path}/{file_prefix}{ts_prefix}.{date_suffix}')
       A.append(ds_ts.timeMonthly_avg_areaCellGlobal.values[0])
       V.append(ds_ts.timeMonthly_avg_volumeCellGlobal.values[0])
       accumulated_ice_runoff_flux.append(
           ds.accumulatedIceRunoffFlux.values[0])
       accumulated_river_runoff_flux.append(
           ds.accumulatedRiverRunoffFlux.values[0])
       accumulated_removed_ice_runoff_flux.append(
           -ds.accumulatedRemovedIceRunoffFlux.values[0])
       accumulated_removed_river_runoff_flux.append(
           -ds.accumulatedRemovedRiverRunoffFlux.values[0])
       accumulated_iceberg_flux.append(
           ds.accumulatedIcebergFlux.values[0])
       accumulated_land_ice_flux.append(
           ds.accumulatedLandIceFlux.values[0])
       accumulated_rain_flux.append(
           ds.accumulatedRainFlux.values[0])
       accumulated_snow_flux.append(
           ds.accumulatedSnowFlux.values[0])
       accumulated_evaporation_flux.append(
           ds.accumulatedEvaporationFlux.values[0])
       accumulated_sea_ice_flux.append(
           ds.accumulatedSeaIceFlux.values[0])
       accumulated_frazil_flux.append(
           ds.accumulatedFrazilFlux.values[0])
A = np.array(A)
V = np.array(V)
total_ssh = V / A

# SSH rates are given in mm/year
factor = mass_flux_to_SLR * s_year / A
ssh_rate_rain = np.array(accumulated_rain_flux) * factor
ssh_rate_snow = np.array(accumulated_snow_flux) * factor
ssh_rate_evaporation = np.array(accumulated_evaporation_flux) * factor
ssh_rate_sea_ice = np.array(accumulated_sea_ice_flux) * factor
ssh_rate_river_runoff = np.array(accumulated_river_runoff_flux) * factor
ssh_rate_ice_runoff = np.array(accumulated_ice_runoff_flux) * factor
ssh_rate_iceberg = np.array(accumulated_iceberg_flux) * factor
ssh_rate_frazil = np.array(accumulated_frazil_flux) * factor
ssh_rate_land_ice = np.array(accumulated_land_ice_flux) * factor
ssh_rate_removed_river_runoff = np.array(accumulated_removed_river_runoff_flux) * factor
ssh_rate_removed_ice_runoff = np.array(accumulated_removed_ice_runoff_flux) * factor

ssh_rate_land_ice_imbalance = ssh_rate_iceberg + \
                              ssh_rate_land_ice + \
                              ssh_rate_removed_river_runoff + \
                              ssh_rate_removed_ice_runoff
ssh_rate_total = ssh_rate_rain + \
                 ssh_rate_snow + \
                 ssh_rate_evaporation + \
                 ssh_rate_river_runoff + \
                 ssh_rate_ice_runoff + \
                 ssh_rate_sea_ice + \
                 ssh_rate_iceberg + \
                 ssh_rate_frazil + \
                 ssh_rate_land_ice

ssh_anomaly = (total_ssh - total_ssh[0]) * 1000. # convert from m to mm
dt = 1 / 12 # in years
ssh_rate_derived = (ssh_anomaly[1:] - ssh_anomaly[:-1]) / dt # [m] [1/s] * [mm/m] * [s/year]
ssh_anomaly_cum = np.zeros((len(ssh_rate_total)))
ssh_anomaly_land_ice_imbalance = np.zeros((len(ssh_rate_total)))
ssh_start_cum = 0
ssh_start_imbalance = 0
for i in range(len(ssh_anomaly_cum)):
    ssh_anomaly_cum[i] = ssh_start_cum + ssh_rate_total[i] / 12
    ssh_anomaly_land_ice_imbalance[i] = ssh_start_imbalance + ssh_rate_land_ice_imbalance[i] / 12
    ssh_start_cum = ssh_anomaly_cum[i]
    ssh_start_imbalance = ssh_anomaly_land_ice_imbalance[i]
ssh_anomaly_cum = ssh_anomaly_cum - ssh_anomaly_cum[0]
ssh_anomaly_land_ice_imbalance = ssh_anomaly_land_ice_imbalance - ssh_anomaly_land_ice_imbalance[0]
mask = ssh_anomaly != 0.
error_factor = np.mean(np.abs((ssh_anomaly_cum[mask] - ssh_anomaly[mask])/ssh_anomaly[mask]))

fig = plt.figure()
plt.plot(decimal_year, accumulated_land_ice_flux, label='land ice', color='navy')
plt.plot(decimal_year, accumulated_iceberg_flux, label='iceberg', color='lightskyblue')
plt.plot(decimal_year, accumulated_removed_river_runoff_flux, label='removed river runoff', color='orangered')
plt.plot(decimal_year, accumulated_removed_ice_runoff_flux, label='removed ice runoff', color='lightsalmon')
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.legend()
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('d(mass)/dt (kg/s)')
fig.savefig(f'{output_dir}/ssh_flux_components.png')

fig = plt.figure()
plt.plot(decimal_year, ssh_rate_land_ice_imbalance, label='land ice imbalance')
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('dSSH/dt (mm/year)')
fig.savefig(f'{output_dir}/ssh_rate_imbalance.png')

fig = plt.figure()
plt.plot(decimal_year[1:], ssh_rate_derived[:], 'k', label='total')
plt.plot(decimal_year, ssh_rate_land_ice_imbalance, label='land ice imbalance')
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
#plt.plot(decimal_year[1:], ssh_rate_derived[:], label='from volume')
#plt.plot(decimal_year, ssh_rate_total, label='from flux')
plt.legend()
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('dSSH/dt (mm/year)')
fig.savefig(f'{output_dir}/ssh_rate_imbalance_total.png')

fig = plt.figure()
plt.plot(decimal_year, ssh_rate_land_ice + ssh_rate_iceberg, label='land ice + iceberg', color='dodgerblue')
plt.plot(decimal_year, ssh_rate_snow + ssh_rate_rain + ssh_rate_evaporation,
         label='precipitation - evaporation', color='lightgreen')
# plt.plot(decimal_year, ssh_rate_snow + ssh_rate_rain, label='precipitation')
# plt.plot(decimal_year, ssh_rate_evaporation, label='evaporation')
plt.plot(decimal_year, ssh_rate_sea_ice + ssh_rate_frazil, label='sea ice + frazil', color='orchid')
plt.plot(decimal_year, ssh_rate_river_runoff + ssh_rate_ice_runoff,
         label='river runoff + ice runoff', color='orange')
plt.plot(decimal_year[1:], ssh_rate_derived[:], 'k', linewidth=0.5, label='total')
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.legend()
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('dSSH/dt (mm/year)')
fig.savefig(f'{output_dir}/ssh_rate_components.png')

fig = plt.figure()
plt.plot(decimal_year, ssh_rate_land_ice, label='land ice', color='navy')
plt.plot(decimal_year, ssh_rate_iceberg, label='iceberg', color='grey')
plt.plot(decimal_year, ssh_rate_removed_river_runoff, label='removed river runoff', color='orangered')
plt.plot(decimal_year, ssh_rate_removed_ice_runoff, label='removed ice runoff', color='lightsalmon')
plt.plot(decimal_year, ssh_rate_land_ice_imbalance, label='land ice imbalance', color='dodgerblue')
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.legend()
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('dSSH/dt (mm/year)')
fig.savefig(f'{output_dir}/ssh_rate_ice_components.png')

fig = plt.figure()
plt.plot(decimal_year, ssh_anomaly, label='volume-derived')
plt.plot(decimal_year, ssh_anomaly_cum, label='flux-derived')
plt.legend()
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('SSH (mm)')
fig.savefig(f'{output_dir}/ssh_time_series.png')

fig = plt.figure()
plt.plot(decimal_year, ssh_anomaly_land_ice_imbalance,
         color='dodgerblue', label='land_ice_imbalance')
plt.fill_between(decimal_year,
                 (1 - error_factor) * ssh_anomaly_land_ice_imbalance,
                 (1 + error_factor) * ssh_anomaly_land_ice_imbalance,
                 color='dodgerblue', edgecolor=None, alpha=0.5)
plt.legend()
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('SSH (mm)')
fig.savefig(f'{output_dir}/ssh_time_series_imbalance.png')

fig = plt.figure()
plt.plot(decimal_year, ssh_anomaly, 'k', label='total')
plt.plot(decimal_year, ssh_anomaly_land_ice_imbalance,
         color='dodgerblue', label='land_ice_imbalance')
plt.fill_between(decimal_year,
                 (1 - error_factor) * ssh_anomaly_land_ice_imbalance,
                 (1 + error_factor) * ssh_anomaly_land_ice_imbalance,
                 color='dodgerblue', edgecolor=None, alpha=0.5)
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.legend()
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('SSH (mm)')
fig.savefig(f'{output_dir}/ssh_time_series_imbalance_total.png')
