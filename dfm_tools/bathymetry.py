import numpy as np
import xarray as xr


def write_bathy_toasc(filename_asc,lon_sel_ext,lat_sel_ext,elev_sel_ext,asc_fmt='%9.2f',nodata_val=32767):
    """
    empty docstring
    """

    print('writing to asc file')
    if elev_sel_ext.shape != (lat_sel_ext.shape[0], lon_sel_ext.shape[0]):
        raise ValueError('resulting shape of elev_sel_ext %s is not consistent with lat_sel_ext/lon_sel_ext vars %s'%(elev_sel_ext.shape,(lat_sel_ext.shape[0], lon_sel_ext.shape[0])))
    if isinstance(elev_sel_ext,np.ma.core.MaskedArray): #masked has to be filled in order for the nans to be visible
        elev_sel_ext = elev_sel_ext.filled(np.nan)
    if np.isnan(elev_sel_ext).sum()>0:
        elev_sel_ext = elev_sel_ext.copy()
        elev_sel_ext[np.isnan(elev_sel_ext)] = nodata_val
        print('replaced nan values with %d'%(nodata_val))
    resinv_lonlat = np.round(1/np.diff(lon_sel_ext[:2]),2)[0]
    resinv_lat = np.round(1/np.diff(lat_sel_ext[:2]),2)[0]
    if resinv_lonlat!=resinv_lat:
        raise ValueError('inconsistent resolution over lat/lon')
    
    with open(filename_asc,'w') as file_asc:
        file_asc.write('ncols         %d\n'%(elev_sel_ext.shape[1]))
        file_asc.write('nrows         %d\n'%(elev_sel_ext.shape[0]))
        file_asc.write('xllcenter     %13.8f\n'%(lon_sel_ext[0]))
        file_asc.write('yllcenter     %13.8f\n'%(lat_sel_ext[0]))
        file_asc.write('cellsize      %16.11f\n'%(1/resinv_lonlat))
        #file_asc.write('NODATA_value  %d\n'%(nodata_val))
        file_asc.write('NODATA_value  '+asc_fmt%(nodata_val)+'\n')
    with open(filename_asc,'a') as file_asc:
        np.savetxt(file_asc,np.flip(elev_sel_ext,axis=0),fmt=asc_fmt)
    print('...finished')


def read_asc(file_asc:str) -> xr.Dataset:
    """
    Reading asc file into a xarray.Dataset

    Parameters
    ----------
    file_asc : str
        asc file with header ncols, nrows, xllcenter, yllcenter, cellsize, NODATA_value.

    Returns
    -------
    ds_asc : xr.Dataset
        xarray.Dataset with the data from the asc file as an array and 
        the lat/lon coordinates as separate coordinate variables.

    """
    # read header
    header_dict = {}
    with open(file_asc) as f:
        for linenum, line in enumerate(f, 1):
            linesplit = line.split()
            linestart = linesplit[0]
            linestop = linesplit[-1]
            
            try:
                # try to convert it to float, this will fail for headerlines
                # it will succeed for numeric lines, but then the loop breaks
                float(linestart)
                break
            except ValueError:
                # convert header values to int if possible or float if not
                try:
                    header_value = int(linestop)
                except ValueError:
                    header_value = float(linestop)
                header_dict[linestart] = header_value
                skiprows = linenum
    
    # read data
    asc_np = np.loadtxt(file_asc, skiprows=skiprows, dtype=float)
    
    # derive x/y values and assert shape
    num_x = header_dict['ncols']
    num_y = header_dict['nrows']
    start_x = header_dict['xllcenter']
    start_y = header_dict['yllcenter']
    step = header_dict['cellsize']
    nodata = header_dict['NODATA_value']
    
    assert asc_np.shape == (num_y, num_x)
    
    x_vals = np.arange(0, num_x) * step + start_x
    y_vals = np.arange(0, num_y) * step + start_y
    
    # flip over latitude and replace nodata with nan
    asc_np = np.flipud(asc_np)
    asc_np[asc_np==nodata] = np.nan
    
    ds_asc = xr.Dataset()
    ds_asc['lon'] = xr.DataArray(x_vals, dims=('lon'))
    ds_asc['lat'] = xr.DataArray(y_vals, dims=('lat'))
    ds_asc['data'] = xr.DataArray(asc_np, dims=('lat','lon'))
    ds_asc = ds_asc.assign_attrs(header_dict)
    return ds_asc

