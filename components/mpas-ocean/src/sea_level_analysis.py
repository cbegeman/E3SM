import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# case_name = 'v3.LR.CRYO1850-DISMF.paolo'
case_name = 'v3.LR.CRYO1850-DISMF'
file_prefix = f'{case_name}.chrysalis'
#file_prefix = '20240201.v3.SORRM.piControl.chrysalis'
run_path = f'/lcrc/group/acme/ac.dcomeau/E3SMv3/{file_prefix}/archive/ocn/hist/'
#run_path = f'/lcrc/group/acme/ac.dcomeau/E3SMv3/{file_prefix}/run/'
#run_path = '/lcrc/group/acme/ac.dcomeau/E3SMv3_dev/20240201.v3.SORRM.piControl.chrysalis/run/'
conserv_prefix = 'mpaso.hist.am.conservationCheck'
ts_prefix = 'mpaso.hist.am.timeSeriesStatsMonthly'
output_dir = f'/lcrc/group/acme/ac.cbegeman/ssh_analysis/{case_name}'

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
# fluxes have units of kg m^-2 s^-1
# accumulated values have units of kg s^-1
for year in range(5):
   for month in range(12):
       decimal_year.append(year + 1. + month / 12.)
       date_suffix = f'{year+1:04g}-{month+1:02g}-01.nc'
       ds = xr.open_dataset(f'{run_path}/{file_prefix}.{conserv_prefix}.{date_suffix}')
       ds_ts = xr.open_dataset(f'{run_path}/{file_prefix}.{ts_prefix}.{date_suffix}')
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

print(accumulated_land_ice_flux)
print(f'dismf: {np.mean(accumulated_land_ice_flux[1:]) * s_year * 1e-12} Gt/yr')
#print(np.mean(accumulated_iceberg_flux + \
#              accumulated_land_ice_flux + \
#              accumulated_removed_river_runoff_flux + \
#              accumulated_removed_ice_runoff_flux) * s_year / np.mean(A))
#print(np.mean(accumulated_rain_flux + \
#              accumulated_snow_flux + \
#              accumulated_evaporation_flux + \
#              accumulated_river_runoff_flux + \
#              accumulated_ice_runoff_flux + \
#              accumulated_sea_ice_flux + \
#              accumulated_iceberg_flux + \
#              accumulated_frazil_flux + \
#              accumulated_land_ice_flux) * s_year / np.mean(A))
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
ssh_rate_data_land_ice = 1.712 + 0.867
ssh_rate_steady_land_ice = 1100/361.8 # same as 1994 rate
ssh_rate_high_land_ice = 1500/361.8
ssh_rate_obs_land_ice = (2660/(2017 - 1992))/361.8
ssh_rate_removed_river_runoff = np.array(accumulated_removed_river_runoff_flux) * factor
ssh_rate_removed_ice_runoff = np.array(accumulated_removed_ice_runoff_flux) * factor
#print(np.mean(ssh_rate_iceberg + ssh_rate_land_ice))
#print(np.mean(ssh_rate_removed_ice_runoff + ssh_rate_removed_river_runoff))

print(f'ssh_rate_land_ice: {np.mean(ssh_rate_land_ice)}')
mass_flux_land_ice_imbalance = (np.array(accumulated_iceberg_flux) + \
                               np.array(accumulated_land_ice_flux) + \
                               np.array(accumulated_removed_river_runoff_flux) + \
                               np.array(accumulated_removed_ice_runoff_flux)) * s_year * 1e-12
mass_anomaly_land_ice_imbalance = np.cumsum(mass_flux_land_ice_imbalance) * (s_year/12)
ssh_rate_land_ice_imbalance = ssh_rate_iceberg + \
                              ssh_rate_land_ice + \
                              ssh_rate_removed_river_runoff + \
                              ssh_rate_removed_ice_runoff
ssh_rate_land_ice_imbalance_corrected = ssh_rate_iceberg + \
                              ssh_rate_data_land_ice * np.ones((len(ssh_rate_iceberg))) + \
                              ssh_rate_removed_river_runoff + \
                              ssh_rate_removed_ice_runoff
print(f'ssh_rate_land_ice_imbalance: {np.mean(ssh_rate_land_ice_imbalance)}')
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

#print(np.mean(ssh_rate_derived) * 1e-3 * rho_sw * np.mean(A)) # [mm/yr] * [m/mm] * [kg/m^3] * [m^2]

ssh_anomaly_cum = np.zeros((len(ssh_rate_total)))
ssh_anomaly_land_ice_imbalance = np.zeros((len(ssh_rate_total)))
ssh_anomaly_land_ice_imbalance_corrected = np.zeros((len(ssh_rate_total)))
ssh_start_cum = 0
ssh_start_imbalance = 0
ssh_start_imbalance_corrected = 0
for i in range(len(ssh_anomaly_cum)):
    ssh_anomaly_cum[i] = ssh_start_cum + ssh_rate_total[i] / 12
    ssh_anomaly_land_ice_imbalance[i] = ssh_start_imbalance + ssh_rate_land_ice_imbalance[i] / 12
    ssh_anomaly_land_ice_imbalance_corrected[i] = ssh_start_imbalance_corrected + ssh_rate_land_ice_imbalance_corrected[i] / 12
    ssh_start_cum = ssh_anomaly_cum[i]
    ssh_start_imbalance = ssh_anomaly_land_ice_imbalance[i]
    ssh_start_imbalance_corrected = ssh_anomaly_land_ice_imbalance_corrected[i]
ssh_anomaly_cum = ssh_anomaly_cum - ssh_anomaly_cum[0]
ssh_anomaly_land_ice_imbalance = ssh_anomaly_land_ice_imbalance - ssh_anomaly_land_ice_imbalance[0]
ssh_anomaly_land_ice_imbalance_corrected = ssh_anomaly_land_ice_imbalance_corrected - ssh_anomaly_land_ice_imbalance_corrected[0]
mask = ssh_anomaly != 0.
#error_factor = np.mean(np.abs((ssh_rate_total[1:] - ssh_rate_derived)/ssh_rate_derived))
#error_factors = np.abs((ssh_anomaly[mask] - ssh_anomaly_cum[mask])/ssh_anomaly[mask])
#error_factor = np.mean(error_factors[100:])
#print(error_factor)
#print(np.min(error_factors),np.max(error_factors))

fig = plt.figure()
plt.plot(decimal_year, accumulated_land_ice_flux, label='land ice', color='navy')
plt.plot(decimal_year, accumulated_iceberg_flux, label='iceberg', color='lightskyblue')
plt.plot(decimal_year, accumulated_removed_river_runoff_flux, label='removed river runoff', color='orangered')
plt.plot(decimal_year, accumulated_removed_ice_runoff_flux, label='removed ice runoff', color='lightsalmon')
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('d(mass)/dt (kg/s)')
fig.savefig(f'{output_dir}/ssh_flux_components.png', bbox_inches="tight")

fig = plt.figure()
plt.plot(decimal_year, mass_flux_land_ice_imbalance, label='land ice imbalance')
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('dM/dt (Gt/yr)')
fig.savefig(f'{output_dir}/mass_flux_imbalance.png', bbox_inches="tight")

fig = plt.figure()
plt.plot(decimal_year, mass_anomaly_land_ice_imbalance, label='land ice imbalance')
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('M (Gt)')
fig.savefig(f'{output_dir}/mass_anomaly_imbalance.png', bbox_inches="tight")

fig = plt.figure()
plt.plot(decimal_year, ssh_rate_land_ice_imbalance, label='land ice imbalance')
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('dSSH/dt (mm/year)')
fig.savefig(f'{output_dir}/ssh_rate_imbalance.png', bbox_inches="tight")

fig = plt.figure()
plt.plot(decimal_year[1:], ssh_rate_derived[:], 'k', label='total')
plt.plot(decimal_year, ssh_rate_land_ice_imbalance, label='land ice imbalance')
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
#plt.plot(decimal_year[1:], ssh_rate_derived[:], label='from volume')
#plt.plot(decimal_year, ssh_rate_total, label='from flux')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('dSSH/dt (mm/year)')
fig.savefig(f'{output_dir}/ssh_rate_imbalance_total.png', bbox_inches="tight")

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
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('dSSH/dt (mm/year)')
fig.savefig(f'{output_dir}/ssh_rate_components.png', bbox_inches="tight")

fig = plt.figure()
plt.plot(decimal_year, ssh_rate_land_ice, label='land ice', color='navy')
plt.plot(decimal_year, ssh_rate_iceberg, label='iceberg', color='grey')
plt.plot(decimal_year,
         ssh_rate_removed_river_runoff + ssh_rate_removed_ice_runoff,
         label='removed river + ice runoff', color='lightsalmon')
#plt.plot(decimal_year, ssh_rate_removed_river_runoff, label='removed river runoff', color='orangered')
#plt.plot(decimal_year, ssh_rate_removed_ice_runoff, label='removed ice runoff', color='lightsalmon')
#plt.plot(decimal_year, ssh_rate_land_ice_imbalance, label='land ice imbalance', color='dodgerblue')
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('dSSH/dt (mm/year)')
fig.savefig(f'{output_dir}/ssh_rate_ice_components_reduced.png', bbox_inches="tight")

fig = plt.figure()
plt.plot(decimal_year, ssh_anomaly, label='volume-derived')
plt.plot(decimal_year, ssh_anomaly_cum, label='flux-derived')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('SSH (mm)')
fig.savefig(f'{output_dir}/ssh_time_series.png', bbox_inches="tight")


fig = plt.figure()
plt.title(case_name)
plt.plot(decimal_year, ssh_anomaly_land_ice_imbalance,
         color='dodgerblue', label='simulated AIS imbalance')
#plt.fill_between(decimal_year,
#                 (1 - error_factor) * ssh_anomaly_land_ice_imbalance,
#                 (1 + error_factor) * ssh_anomaly_land_ice_imbalance,
#                 color='dodgerblue', edgecolor=None, alpha=0.5)
#plt.fill_between(decimal_year,
#                 (1 - error_factor) * ssh_anomaly_land_ice_imbalance,
#                 ssh_anomaly_land_ice_imbalance,
#                 color='dodgerblue', edgecolor=None, alpha=0.5)
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('SSH anomaly (mm)')
fig.savefig(f'{output_dir}/ssh_time_series_imbalance.png', bbox_inches="tight")

plt.plot(decimal_year, ssh_anomaly,
         'k', label='simulated sea level anomaly', linewidth=0.5)
fig.savefig(f'{output_dir}/ssh_time_series_imbalance_cum.png', bbox_inches="tight")

plt.plot(decimal_year,
         ssh_anomaly_land_ice_imbalance_corrected,
         color='lightskyblue', label='corrected AIS imbalance')
#plt.plot([min(decimal_year), max(decimal_year)],
#         [0, ssh_rate_data_land_ice * max(decimal_year)],
#         color='orange', label='data land ice melt')
plt.plot([min(decimal_year), max(decimal_year)],
         [0, ssh_rate_obs_land_ice * max(decimal_year)],
         color='green', label='observed AIS imbalance')
#plt.fill_between([min(decimal_year), max(decimal_year)],
#                 [0, ssh_rate_steady_land_ice * max(decimal_year)],
#                 [0, ssh_rate_high_land_ice * max(decimal_year)],
#                 label='observed land ice melt',
#                 color='lightgreen', edgecolor=None, alpha=0.5)
plt.plot([min(decimal_year), max(decimal_year)], [0,0], ':k')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlim([min(decimal_year), max(decimal_year)])
plt.xlabel('Year')
plt.ylabel('SSH (mm)')
fig.savefig(f'{output_dir}/ssh_time_series_imbalance_corrected.png', bbox_inches="tight")

plt.plot(decimal_year, ssh_anomaly_cum, color='lightsalmon',
         label='flux-derived sea level anomaly', linewidth=0.5)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig.savefig(f'{output_dir}/ssh_time_series_imbalance_corrected_cum.png', bbox_inches="tight")

#plt.fill_between(decimal_year,
#                 (1 - error_factor) * ssh_anomaly_land_ice_imbalance,
#                 ssh_anomaly_land_ice_imbalance,
#                 color='dodgerblue', edgecolor=None, alpha=0.5)
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#fig.savefig(f'{output_dir}/ssh_time_series_imbalance_corrected_uncertainty.png', bbox_inches="tight")
