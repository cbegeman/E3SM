import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

s_year = 3.154e7 
rho_sw = 1026.
mass_flux_to_SLR = 1.e3 / rho_sw # we still need to divide by area of the ocean below to get units of mm/s

#dismf_prefix = 'prescribed_ismf_adusumilli2020'
dismf_prefix = 'prescribed_ismf_paolo2023'

dismf_area_cell = 2.5e5  # m^2 because the dataset is on a 500m grid
# The MPAS-Analysis source file does not report cell areas
dismf_ref = '/lcrc/group/e3sm/public_html/diagnostics/observations/Ocean/Melt/Adusumilli/Adusumilli_2020_iceshelf_melt_rates_2010-2018_v0.20230504.nc'
#dismf_new = '/lcrc/group/e3sm/data/inputdata/ocn/mpas-o/RRSwISC6to18E3r4/prescribed_ismf_paolo2023.RRSwISC6to18E3r4.20240221.nc'
ds_ref = xr.open_dataset(dismf_ref)
ref_melt_rate = ds_ref['meltRate'].values # m/yr
total_melt_volume_flux = np.nansum(ref_melt_rate * dismf_area_cell) # m^3/yr
total_melt_mass_flux = total_melt_volume_flux * rho_sw # kg/yr
# print(f'Adusumilli: {np.sum(total_melt_mass_flux) * 1e-15} x 10^15 kg/yr')
# This file is based on 2010-2018 observations
print(f'Adusumilli: {np.sum(total_melt_mass_flux) * 1e-12} Gt/yr')
# print('Adusumilli: 1260 ± 150 Gt/yr')
# Adusumilli et al. 2020: The total meltwater flux, based on the area-
# integrated basal melt rate over all Antarctic ice shelves averaged over
# 1994–2018, was 1,260 ± 150 Gt yr–1

base_path = '/lcrc/group/e3sm/data/inputdata/ocn/mpas-o'

mesh_dict = {'HR': {'mesh_name': 'RRSwISC6to18E3r4',
                    'dismf_suffix': '20240221',
                    'mesh_suffix': '20240105'},
             'SORRM': {'mesh_name': 'SOwISC12to30E3r2',
                       'dismf_suffix': '20240221',
                       'mesh_suffix': '20240104'}
            }

for entry in mesh_dict.keys():
    mesh_name = mesh_dict[entry]['mesh_name']
    mesh_suffix = mesh_dict[entry]['mesh_suffix']
    dismf_suffix = mesh_dict[entry]['dismf_suffix']
    dismf_file = f'{base_path}/{mesh_name}/{dismf_prefix}.{mesh_name}.{dismf_suffix}.nc'
    mesh_file = f'{base_path}/{mesh_name}/mpaso.{mesh_name}.{mesh_suffix}.nc'
    ds_mesh = xr.open_dataset(mesh_file)
    ds = xr.open_dataset(dismf_file)

    # fluxes have units of kg m^-2 s^-1
    area_cell = ds_mesh['areaCell'].values
    valid_cell = ds_mesh['maxLevelCell'].values > 0
    data_fw_flux = ds['dataLandIceFreshwaterFlux'].values

    # First we just convert to kg yr^-1
    data_fw_flux_total = data_fw_flux * area_cell * s_year
    
    data_fw_SLR = np.sum(data_fw_flux_total) * mass_flux_to_SLR / np.sum(area_cell[valid_cell])
    
    # print(f'{entry}: {np.sum(data_fw_flux_total) * 1e-15} x 10^15 kg/yr')
    print(f'{entry}: {np.sum(data_fw_flux_total) * 1e-12} Gt/yr')
    print(f'{entry}: {data_fw_SLR} mm/yr')
