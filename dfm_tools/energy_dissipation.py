import datetime as dt
from dask.diagnostics import ProgressBar

__all__ = ["compute_energy_dissipation"]

def compute_energy_dissipation(data_xr_map,file_ED_computed):
    """
    Example:
        data_xr_map = dfmt.open_partitioned_dataset(file_nc_map,chunks={'time':100}) #TODO: important to have time>1, otherwise time-mean floods memory (100-200 seems optimal for GTSM 1month, but memory still floods so 1 year would probably be impossible).
        data_xr_map = data_xr_map.sel(time=slice('2014-01-01','2014-02-01'))
        dfmt.compute_energy_dissipation(data_xr_map,file_ED_computed)

    
    TODO: clean up these comments
    
    Stel de bodemwrijving is (tau_x,tau_y) dan is het verlies aan energie in W/m^2
    E_loss = tau_x*u_x + tau_y*u_y
    
    Voor Chezy geldt (zo uit mijn hoofd):
    tau_x = - (rho*g)/(C^2) u_x * sqrt(u_x^2+u_y^2)
    tau_y = - (rho*g)/(C^2) u_y * sqrt(u_x^2+u_y^2)
    
    En gecombineerd
    E_loss = - (rho*g)/(C^2) * (u_x^2+u_y^2)^(3/2)
    
    Approximation: using velocities in cell centers, this is averaged but easiest
    
    E_loss(area) = sum_cells [ - (rho*g)/(C^2) * (sqrt(u_x^2+u_y^2))^3 ] * cell_area
    
    Dit is dan voor een tijdstip. Ik zou middelen over een iets langere periode zodat het de variatie over het getij verdwijnt, want de waarde is bij springtij waarschijnlijk groter.

    """
    rho = 1020 #kg/m3
    g = 9.81 #m/s2
    
    attrs_ED = {'long_name':'energy dissipation','units':'W'}    
    attrs_ED_pm2 = {'long_name':'energy dissipation','units':'W/m^2'}    
    attrs_ED_areasum = {'long_name':'area-sum energy dissipation', 'units':'W'}
    attrs_ED_timemean = {'long_name':'springneaptide-mean energy dissipation', 'units': 'W'}
    attrs_ED_pm2_timemean = {'long_name':'springneaptide-mean energy dissipation', 'units':'Wm^2'}
    
    print('>> compute energyloss from C/u/area: ',end='')
    dtstart = dt.datetime.now()
    data_ucmag = data_xr_map.mesh2d_ucmag #m/s
    data_czs = data_xr_map.mesh2d_czs #m0.5s-1
    data_czs = data_czs.where(data_czs!=0) #0 will result in inf energyloss, so replace by nan
    data_ba = data_xr_map.mesh2d_flowelem_ba #m2
    data_xr_map['ED'] = ((rho*g)/(data_czs**2) * data_ucmag**3 * data_ba).assign_attrs(attrs_ED) #equation: (kg/m3)*(m/s2)/(m/s2)*(m3/s3)*m2 = kg/s3*m2 = W = J/s
    data_xr_map['ED_pm2'] = (data_xr_map['ED'] / data_ba).assign_attrs(attrs_ED_pm2) #TODO: maybe include ba instead of supplying _pm2 separately?
    
    #TODO: check da.sum(split_every=4): https://github.com/dask/dask/issues/883. Does not seem to work: "TypeError: nansum() got an unexpected keyword argument 'split_every'"
    #TODO: chunking: float32 = 32 bits = 4 bytes >> use this to calculate chunking sizes (100-200MB per chunk is optimal)
    data_xr_map['ED_areasum'] = data_xr_map.ED.sum(dim='mesh2d_nFaces').assign_attrs(attrs_ED_areasum)
    data_xr_map['ED_timemean'] = data_xr_map.ED.mean(dim='time').assign_attrs(attrs_ED_timemean)
    data_xr_map['ED_pm2_timemean'] = data_xr_map.ED_pm2.mean(dim='time').assign_attrs(attrs_ED_pm2_timemean)
    
    
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    data_xr_map_computed = data_xr_map[['mesh2d_flowelem_ba','ED_areasum','ED_timemean','ED_pm2_timemean']] #has to contain some var with cells (not only ED_areasum), otherwise ugrid accessor is not valid
    #data_xr_map_computed.ED_timemean.data.visualize() #does not yet work >> visualizes tasks of dask array
    
    print('>> compute and save data_xr_map_computed to netcdf: ')
    with ProgressBar():
        data_xr_map_computed.ugrid.to_netcdf(file_ED_computed)
    
    data_xr_map_computed.close()
    