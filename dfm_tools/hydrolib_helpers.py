import pandas as pd
import cftime
import numpy as np
import xarray as xr
from cftime import date2num
import hydrolib.core.dflowfm as hcdfm
import datetime as dt
import geopandas
from shapely.geometry import LineString

__all__ = ["Dataset_to_T3D",
           "Dataset_to_TimeSeries",
           "Dataset_to_Astronomic",
           "DataFrame_to_PolyObject",
           "geodataframe_to_PolyFile",
           "DataFrame_to_TimModel",
           "ForcingModel_to_plipointsDataset",
           "forcinglike_to_Dataset",
           "pointlike_to_DataFrame",
           "TimModel_to_DataFrame",
           "pointlike_to_geodataframe_points",
           "PolyFile_to_geodataframe_points",
           "PolyFile_to_geodataframe_linestrings",
           "gdf_linestrings_to_points",
           "tekalobject_to_DataFrame",
       ]


def get_ncbnd_construct():
    """
    netcdf structure is based on FEWS example file:
    p:\\dflowfm\\maintenance\\JIRA\\06000-06999\\06187\\C01_JV\\salinity_DCSM-FM_OB_all.nc
    Some improvement suggestions are provided in:
    p:\\dflowfm\\maintenance\\JIRA\\06000-06999\\06187\\C01_JV\\convert_netcdf_bnd.py
    """
    attrs_pointx = {"standard_name": "longitude",
                    "long_name": "longitude",
                    "units": "degrees_east",
                    "axis": "X"
                    }
    attrs_pointy = {"standard_name": "latitude",
                    "long_name": "latitude",
                    "units": "degrees_north",
                    "axis": "Y"
                    }
    attrs_depth = {"standard_name":"z",
                   "long_name":"z",
                   "units":"m",
                   "positive":"up",
                   "axis":"Z",
                    }
    
    ncbnd_construct = {"varn_depth":"z",
                       "dimn_depth":"z",
                       "varn_pointx":"lon",
                       "varn_pointy":"lat",
                       "varn_pointname":"station_id",
                       "dimn_point":"node",
                       "attrs_pointx":attrs_pointx,
                       "attrs_pointy":attrs_pointy,
                       "attrs_depth":attrs_depth,
                       }
    
    return ncbnd_construct


def Dataset_to_T3D(datablock_xr):
    """
    convert an xarray.DataArray (is one data_var) or an xarray.Dataset (with one or two data_vars) with time and depth dimension to a hydrolib T3D object
    """
    
    ncbnd_construct = get_ncbnd_construct()
    dimn_depth = ncbnd_construct['dimn_depth']
    varn_depth = ncbnd_construct['varn_depth']
    
    if not isinstance(datablock_xr,(xr.DataArray,xr.Dataset)):
        raise TypeError(f'expected xarray.DataArray or xarray.Dataset, not {type(datablock_xr)}')
        
    vector = False
    if isinstance(datablock_xr,xr.DataArray):
        data_xr_var0 = datablock_xr
    elif isinstance(datablock_xr,xr.Dataset):
        data_vars = list(datablock_xr.data_vars)
        data_xr_var0 = datablock_xr[data_vars[0]]
        if len(data_vars)==2:
            if not pd.Series(data_vars).isin(['ux','uy']).all():
                raise Exception(f'Dataset with 2 data_vars should contain only ux/uy data_vars, but contains {data_vars}')
            vector = True
            data_xr_var1 = datablock_xr[data_vars[1]]
        elif len(data_vars) > 2:
            raise ValueError(f'Dataset should contain 1 or 2 data_vars, but contains {len(data_vars)} variables')
    
    #ffill/bfill nan data along over depth dimension (corresponds to vertical extrapolation)
    data_xr_var0 = data_xr_var0.bfill(dim=dimn_depth).ffill(dim=dimn_depth)
    if vector:
        data_xr_var1 = data_xr_var1.bfill(dim=dimn_depth).ffill(dim=dimn_depth)
    
    #TODO: clean up these first lines of code and add description to docstring?
    locationname = data_xr_var0.attrs['locationname']
    refdate_str = data_xr_var0.time.encoding['units']
    
    if not set(data_xr_var0.dims).issubset(set(('time',dimn_depth))): #check if both time and depth dimensions are present
        raise ValueError(f"data_var in provided data_xr has dimensions {data_xr_var0.dims} while ('time',{dimn_depth}) is expected")
    
    #get depth variable and values
    depth_array = data_xr_var0[varn_depth].to_numpy()
    
    #get datablock and concatenate with relative time data
    if vector:
        data_xr_var0_np = data_xr_var0.to_numpy()
        data_xr_var1_np = data_xr_var1.to_numpy()
        datablock_np = np.stack((data_xr_var0_np,data_xr_var1_np),2).reshape(data_xr_var0_np.shape[0],-1) #merge data with alternating rows
    else:
        datablock_np = data_xr_var0.to_numpy()
    
    timevar_sel_rel = date2num(pd.DatetimeIndex(data_xr_var0.time.to_numpy()).to_pydatetime(),units=refdate_str,calendar='standard')
    datablock_incltime = np.concatenate([timevar_sel_rel[:,np.newaxis],datablock_np],axis=1)
    
    # Each .bc file can contain 1 or more timeseries, in this case one for each support point
    verticalpositions_idx = np.arange(data_xr_var0[varn_depth].size)+1
    if vector: #vector T3D object
        QUP_quan_list = [hcdfm.QuantityUnitPair(quantity=quan, unit=data_xr_var0.attrs['units'], vertpositionindex=iVP) for iVP in verticalpositions_idx for quan in data_vars]
        QUP_quan_vector = hcdfm.VectorQuantityUnitPairs(vectorname='uxuyadvectionvelocitybnd', #TODO: vectorname from global attr? (then also support other vectors which is not necessary)
                                                  elementname=data_vars,
                                                  quantityunitpair=QUP_quan_list)
        quantityunitpair = [hcdfm.QuantityUnitPair(quantity="time", unit=refdate_str)]+[QUP_quan_vector]
    else: #normal T3D object
        QUP_quan_list = [hcdfm.QuantityUnitPair(quantity=data_xr_var0.name, unit=data_xr_var0.attrs['units'], vertpositionindex=iVP) for iVP in verticalpositions_idx]
        quantityunitpair=[hcdfm.QuantityUnitPair(quantity="time", unit=refdate_str)]+QUP_quan_list
    
    T3D_object = hcdfm.T3D(name=locationname,
                           #offset=0,
                           #factor=1,
                           vertpositions=depth_array.tolist(),
                           vertinterpolation='linear', #TODO: make these parameters user defined (via attrs)
                           vertPositionType='ZDatum',
                           quantityunitpair=quantityunitpair,
                           timeinterpolation='linear',
                           datablock=datablock_incltime.tolist(),
                           )
    
    return T3D_object


def Dataset_to_TimeSeries(datablock_xr):
    """
    convert an xarray.DataArray or xarray.Dataset with time dimension to a hydrolib TimeSeries object
    """
    
    if not isinstance(datablock_xr,(xr.DataArray,xr.Dataset)):
        raise TypeError(f'Dataset_to_TimeSeries expects xr.DataArray or xr.Dataset, not {type(datablock_xr)}')
    
    if isinstance(datablock_xr,xr.Dataset): #convert Dataset to DataArray
        data_vars = list(datablock_xr.data_vars)
        if len(data_vars)!=1:
            raise ValueError('more than one variable supplied in Dataset, not yet possible') #TODO: add support for multiple quantities and for vectors
        datablock_xr = datablock_xr[data_vars[0]]
    
    #TODO: clean up these first lines of code and add description to docstring?
    locationname = datablock_xr.attrs['locationname']
    bcvarname = datablock_xr.name
    refdate_str = datablock_xr.time.encoding['units']
    
    if datablock_xr.dims != ('time',):
        raise ValueError(f"datablock_xr provided to DataArray_to_TimeSeries has dimensions {datablock_xr.dims} while ('time') is expected")
    
    #get datablock and concatenate with relative time data
    datablock_np = datablock_xr.to_numpy()[:,np.newaxis]
    timevar_sel_rel = date2num(pd.DatetimeIndex(datablock_xr.time.to_numpy()).to_pydatetime(),units=refdate_str,calendar='standard')
    datablock_incltime = np.concatenate([timevar_sel_rel[:,np.newaxis],datablock_np],axis=1)
    
    # Each .bc file can contain 1 or more timeseries, in this case one for each support point
    TimeSeries_object = hcdfm.TimeSeries(name=locationname,
                                         quantityunitpair=[hcdfm.QuantityUnitPair(quantity="time", unit=refdate_str),
                                                           hcdfm.QuantityUnitPair(quantity=bcvarname, unit=datablock_xr.attrs['units'])],
                                         timeinterpolation='linear', #TODO: make userdefined via attrs?
                                         datablock=datablock_incltime.tolist(), 
                                         )
    return TimeSeries_object


def Dataset_to_Astronomic(datablock_xr):
    """
    convert an xarray.Dataset (with amplitude and phase data_vars) to a hydrolib Astronomic object
    
    """
    if not isinstance(datablock_xr,xr.Dataset):
        raise TypeError(f'Dataset_to_Astronomic expects xr.Dataset, not {type(datablock_xr)}')
    
    data_vars = list(datablock_xr.data_vars)
    if 'amplitude' not in data_vars or 'phase_new' not in data_vars:
        raise KeyError('amplitude and/or phase_new not in input xr.Dataset')

    #TODO: clean up these first lines of code and add description to docstring?
    locationname = datablock_xr['amplitude'].attrs['locationname']
        
    #get datablock and concatenate with component names
    datablock_np_cna = datablock_xr['compno'].to_numpy()[:,np.newaxis]
    datablock_np_amp = datablock_xr['amplitude'].to_numpy()[:,np.newaxis]
    datablock_np_phs = datablock_xr['phase_new'].to_numpy()[:,np.newaxis]
    datablock_inclcomp = np.concatenate([datablock_np_cna,datablock_np_amp,datablock_np_phs],axis=1)
    
    Astronomic_object = hcdfm.Astronomic(name=locationname,
                                         quantityunitpair=[hcdfm.QuantityUnitPair(quantity="astronomic component", unit='-'),
                                                           hcdfm.QuantityUnitPair(quantity='waterlevelbnd amplitude', unit=datablock_xr['amplitude'].attrs['units']),
                                                           hcdfm.QuantityUnitPair(quantity='waterlevelbnd phase', unit=datablock_xr['phase'].attrs['units'])],
                                         datablock=datablock_inclcomp.tolist(), 
                                         )
    return Astronomic_object


def DataFrame_to_PolyObject(poly_pd,name,content=None):
    """
    convert a pandas dataframe with x/y columns (and optional others like z/data/comment) to a hydrolib PolyObject
    """
    if 'z' in poly_pd.columns:
        nondata_cols = ['x','y','z']
    else:
        nondata_cols = ['x','y']
    poly_pd_xy = poly_pd[nondata_cols]
    poly_pd_data = pd.DataFrame({'data':poly_pd.drop(nondata_cols,axis=1).values.tolist()})
    poly_pd_polyobj = pd.concat([poly_pd_xy,poly_pd_data],axis=1)
    pointsobj_list = poly_pd_polyobj.T.apply(dict).tolist() #TODO: maybe faster with list iteration
    polyobject = hcdfm.PolyObject(metadata={'name':name,'n_rows':poly_pd.shape[0],'n_columns':poly_pd.shape[1]}, points=pointsobj_list)
    if content is not None:
        polyobject.description = {'content':content}
    return polyobject


def validate_polyline_names(polyfile_obj):
    # TODO: not allowed to have empty or duplicated polyline names in a polyfile, this is not 
    # catched by hydrolib-core: https://github.com/Deltares/HYDROLIB-core/issues/483
    # therefore, we check it here
    names = [x.metadata.name for x in polyfile_obj.objects]
    if len(set(names)) != len(names):
        raise ValueError(f'duplicate polyline names found in polyfile: {names}')
    first_alpha = [x[0].isalpha() for x in names]
    if not all(first_alpha):
        raise ValueError(f'names in polyfile do not all start with a letter: {names}')


def geodataframe_to_PolyFile(poly_gdf, name="L"):
    """
    convert a geopandas geodataframe with x/y columns (and optional others like z/data/comment) to a hydrolib PolyFile
    """
    
    # catch some invalid occurences of name
    if not isinstance(name, str):
        raise TypeError("name should be a string")
    
    # add name column if not present
    if 'name' not in poly_gdf.columns:
        # make a copy to avoid alternating the geodataframe
        poly_gdf = poly_gdf.copy()
        name_nums = poly_gdf.reset_index().index+1
        poly_gdf['name'] = name + name_nums.astype(str)
    
    polyfile_obj = hcdfm.PolyFile()
    # TODO: now only name+geometry, still add other data columns
    for irow, gdf_row in poly_gdf.iterrows():
        poly_geom = gdf_row.geometry
        name_str = gdf_row['name']
        if isinstance(poly_geom, LineString):
            poly_geom_np = np.array(poly_geom.xy).T
        else: # isinstance(poly_geom, shapely.Polygon):
            poly_geom_np = np.array(poly_geom.boundary.xy).T
        pd_columns = ['x','y']
        if poly_geom_np.shape[1]==3:
            pd_columns = ['x','y','z']
        poly_geom_df = pd.DataFrame(poly_geom_np,columns=pd_columns)
        pointsobj_list = poly_geom_df.T.apply(dict).tolist()
        for pnt in pointsobj_list:
            pnt['data'] = []
        polyobject = hcdfm.PolyObject(metadata={'name':name_str,
                                                'n_rows':poly_geom_np.shape[0],
                                                'n_columns':poly_geom_np.shape[1]}, 
                                      points=pointsobj_list)
        polyfile_obj.objects.append(polyobject)
    
    validate_polyline_names(polyfile_obj)
    return polyfile_obj


def DataFrame_to_TimModel(tim_pd, refdate:(dt.datetime, pd.Timestamp, str)):
    """
    converts data from tim_pd to TimModel and puts all headers as comments. Ignores the index, assumes first column is time in minutes since a refdate.
    """
    #TODO: add conversion from datetimes in index to minutes, maybe drop minutes column upon reading? First await https://github.com/Deltares/HYDROLIB-core/issues/511
    
    refdate_pd = pd.Timestamp(refdate)
    
    data_tim = tim_pd.values.tolist()
    times_tim = ((tim_pd.index - refdate_pd).total_seconds()/60).tolist()
    dict_tim = [hcdfm.TimRecord(time=t,data=d) for t,d in zip(times_tim,data_tim)]
    
    comments_datacols = tim_pd.columns.tolist()
    comments = [tim_pd.index.name] + comments_datacols
    for iC, comment in enumerate(comments):
        comments[iC] = f'COLUMN {iC+1}: {comment}'
    timmodel = hcdfm.TimModel(timeseries=dict_tim, comments=comments)
    
    return timmodel


def ForcingModel_to_plipointsDataset(forcingmodel:hcdfm.ForcingModel, npoints=None, convertnan=False) -> xr.Dataset:
    if not isinstance(forcingmodel, hcdfm.ForcingModel):
        raise TypeError('ForcingModel_to_plipointsDataset expects type hcdfm.ForcingModel, not type {type(forcingobj)}')

    ncbnd_construct = get_ncbnd_construct()
    dimn_point = ncbnd_construct['dimn_point']
    varn_pointname = ncbnd_construct['varn_pointname']
    plipointsDataset_list = []
    for forcinglike in forcingmodel.forcing[:npoints]:
        ds_onepoint = forcinglike_to_Dataset(forcinglike, convertnan=convertnan)
        
        #set longname attr
        for datavar in ds_onepoint.data_vars:
            longname = datavar
            if datavar in ['ux','uy']: #TODO: hardcoded behaviour is consitent with maybe_convert_fews_to_dfmt() and elsewhere in dfm_tools, but not desireable
                longname = 'uxuyadvectionvelocitybnd'
            ds_onepoint[datavar] = ds_onepoint[datavar].assign_attrs({'long_name': longname})
        
        # expand pointdim and add pointname as var
        ds_onepoint = ds_onepoint.expand_dims(dimn_point)
        datavar0 = list(ds_onepoint.data_vars)[0]
        pointname = ds_onepoint[datavar0].attrs['locationname']
        ds_onepoint[varn_pointname] = xr.DataArray([pointname],dims=dimn_point)
        ds_onepoint = ds_onepoint.set_coords(varn_pointname)
        
        plipointsDataset_list.append(ds_onepoint)
    ds = xr.concat(plipointsDataset_list, dim=dimn_point)
    
    ds = maybe_convert_fews_to_dfmt(ds)
    
    return ds


def forcinglike_to_Dataset(forcingobj, convertnan=False):
    """
    convert a hydrolib forcing like object (like Timeseries, T3D, Harmonic, etc) to an xarray Dataset with one or more variables.
    
    convertnan: convert depths with the same values over time as the deepest layer to nan (these were created with .bfill() or .ffill()).
    """
    
    ncbnd_construct = get_ncbnd_construct()
    dimn_depth = ncbnd_construct['dimn_depth']
    varn_depth = ncbnd_construct['varn_depth']
    attrs_depth = ncbnd_construct['attrs_depth']
    
    #check if forcingmodel instead of T3D/TimeSeries is provided
    if isinstance(forcingobj, hcdfm.ForcingModel):
        raise TypeError('instead of supplying a ForcingModel, provide a ForcingObject (Timeseries/T3D etc), by doing something like ForcingModel.forcing[0], or use dfmt.ForcingModel_to_plipointsDataset() instead')
    
    allowed_instances = (hcdfm.T3D, hcdfm.TimeSeries, hcdfm.Astronomic)
    if not isinstance(forcingobj, allowed_instances):
        raise TypeError(f'supplied input is not one of: {allowed_instances}')
    
    if isinstance(forcingobj, hcdfm.Astronomic):
        var_quantity_list = [x.quantity for x in forcingobj.quantityunitpair[1:]]
        var_unit = [x.unit for x in forcingobj.quantityunitpair[1:]]
    elif hasattr(forcingobj.quantityunitpair[1],'elementname'): #T3D with vector quantity
        var_quantity_list = forcingobj.quantityunitpair[1].elementname
        var_unit_one = forcingobj.quantityunitpair[1].quantityunitpair[0].unit
        var_unit = [var_unit_one,var_unit_one]
    else: #non-vector TimeSeries or T3D
        var_quantity_list = [forcingobj.quantityunitpair[1].quantity]
        var_unit = [forcingobj.quantityunitpair[1].unit]
    nquan = len(var_quantity_list)
    
    if isinstance(forcingobj, hcdfm.T3D):
        dims = ('time',dimn_depth)
    elif isinstance(forcingobj, hcdfm.TimeSeries):
        dims = ('time')
    elif isinstance(forcingobj, hcdfm.Astronomic):
        dims = ('astronomic_component')
    
    datablock_all = np.array(forcingobj.datablock)
    datablock_data = datablock_all[:,1:] #select all columns except first one (which is the time column)
    if isinstance(forcingobj, hcdfm.Astronomic):
        datablock_data = datablock_data.astype(float) #convert str to float in case of "astronomic component"
    
    data_xr = xr.Dataset()
    for iQ, var_quantity in enumerate(var_quantity_list):
        datablock_data_onequan = datablock_data[:,iQ::nquan]
        datablock_data_onequan = datablock_data_onequan.squeeze() #drop dimensions of len 1 in case of 1 dimension, eg "waterlevelbnd" (first subsetting over depth dimension)
        
        data_xr_var = xr.DataArray(datablock_data_onequan, name=var_quantity, dims=dims)
        if dimn_depth in dims:
            data_xr_var[varn_depth] = forcingobj.vertpositions
            data_xr_var[varn_depth] = data_xr_var[varn_depth].assign_attrs(attrs_depth)
            if convertnan: #convert ffilled/bfilled values back to nan
                deepestlayeridx = data_xr_var[varn_depth].to_numpy().argmin()
                if deepestlayeridx==0: #sorted from deep to shallow layers
                    bool_eq_above = data_xr_var==data_xr_var.shift({varn_depth:-1})
                else: #sorted from shallow to deep layers
                    bool_eq_above = data_xr_var==data_xr_var.shift({varn_depth:1})
                bool_eq_bottom = data_xr_var==data_xr_var.isel({varn_depth:deepestlayeridx})
                bool_nandepths = (bool_eq_above & bool_eq_bottom).all(dim='time')
                data_xr_var = data_xr_var.where(~bool_nandepths)
        if 'time' in dims:
            time_unit = forcingobj.quantityunitpair[0].unit.lower()
            data_xr_var['time'] = cftime.num2pydate(datablock_all[:,0], units=time_unit)
            data_xr_var['time'].encoding['units'] = time_unit #check tz conversion if eg '+01:00' is present in time_unit
            data_xr_var['time'].encoding['calendar'] = 'standard'
        if 'astronomic_component' in dims:
            data_xr_var['astronomic_component'] = datablock_all[:,0]
        
        #add attributes
        data_xr_var.attrs['source'] = 'hydrolib-core object converted to xarray.Dataset with dfm_tools.xarray_helpers.forcinglike_to_Dataset()'
        data_xr_var.attrs['locationname'] = forcingobj.name
        data_xr_var.attrs['units'] = var_unit[iQ]
        forcingobj_keys = forcingobj.__dict__.keys()
        for key in forcingobj_keys: #['comments','name','function','offset','factor','vertinterpolation','vertpositiontype','timeinterpolation']: 
            if key in ['datablock','quantityunitpair','vertpositions','name']: #skipping these since they are in the DataArray already
                continue
            data_xr_var.attrs[key] = str(forcingobj.__dict__[key])
        
        #add DataArray to Dataset
        data_xr[var_quantity] = data_xr_var
    
    return data_xr


def pointlike_to_DataFrame(pointlike, drop_emptycols=True):
    """
    convert a hydrolib object with points (like PolyObject, XYZModel and possibly others) to a pandas DataFrame.

    Parameters
    ----------
    pointlike : TYPE
        Hydrolib-core object with point objects.

    Returns
    -------
    pointlike_pd : TYPE
        DESCRIPTION.
    drop_emptycols : bool, optional
        Drop empty (all-nan) columns automatically, like the z-column in ldb file. The default is True.
    
    Example:
        polyfile_object = PolyFile(file_pli)
        data_pol_pd_list = [pointlike_to_DataFrame(polyobj) for polyobj in polyfile_object.objects]

    """
    
    pointlike_pd = pd.DataFrame([dict(p) for p in pointlike.points])
    if 'data' in pointlike_pd.columns:
        datavals_pd = pd.DataFrame([p.data for p in pointlike.points])
        pointlike_pd = pd.concat([pointlike_pd.drop(['data'],axis=1), datavals_pd],axis=1)
        
    if drop_emptycols:
        pointlike_pd = pointlike_pd.dropna(axis=1).copy()
        
    return pointlike_pd


def TimModel_to_DataFrame(data_tim:hcdfm.TimModel, parse_column_labels:bool = True, refdate:(dt.datetime, pd.Timestamp, str) = None):
    """
    

    Parameters
    ----------
    data_tim : hcdfm.TimModel
        DESCRIPTION.
    parse_column_labels : bool, optional
        Parse column labels from comments. This might fail since there is no standard way of prescribing the columns in the comments. The default is True.
    refdate : (dt.datetime, pd.Timestamp, str, int, float), optional
        DESCRIPTION. The default is None.

    Returns
    -------
    tim_pd : TYPE
        DESCRIPTION.

    """
    #convert to pandas dataframe
    timevals_pd = pd.Index([p.time for p in data_tim.timeseries])
    tim_pd = pd.DataFrame([p.data for p in data_tim.timeseries],index=timevals_pd)
    tim_pd.columns += 2 #make column numbers 1-based, but first column is already in index so start with 2
    tim_pd.index.name = 'time in minutes'
    
    if parse_column_labels:
        #replace column labels with the ones in comments
        tim_pd_columns = tim_pd.columns.tolist()
        for line in data_tim.comments:
            if 'column' in line.lower() and ':' in line: #assume ":" is separator. Remove casing to be able to check for Column/COLUMN/column in string
                sep = ':'
            elif 'column' in line.lower() and '=' in line: #assume "=" is separator. Remove casing to be able to check for Column/COLUMN/column in string
                sep = '='
            else:
                continue
            line_split = line.split(sep)
            colnum = line_split[0].lower().replace('column','').strip()
            if colnum.isnumeric():
                comment_str = ':'.join(line_split[1:]).strip()
                if colnum==1: #time column is now in index, overwrite index name
                    tim_pd.index.name = comment_str
                else:
                    tim_pd_columns[int(colnum)-2] = comment_str
        tim_pd.columns = tim_pd_columns
    
    if refdate:
        refdate_pd = pd.Timestamp(refdate)
        tim_pd.index = refdate_pd + pd.to_timedelta(tim_pd.index,unit='minutes')
    
    return tim_pd


def pointlike_to_geodataframe_points(polyline_object, crs:str=None, add_pointnames=True, only_xy=False):
    """
    empty docstring
    """
    
    ncbnd_construct = get_ncbnd_construct()
    varn_pointname = ncbnd_construct['varn_pointname']
    
    #conversion to dataframe
    polyobject_pd = pd.DataFrame([dict(p) for p in polyline_object.points])
    
    if only_xy:
        df = pd.DataFrame()
    else:
        df = polyobject_pd.drop(['x','y'],axis=1) #also includes n/z and maybe data column
    
    if add_pointnames and hasattr(polyline_object,'metadata'): #optionally add names
        df[varn_pointname] = pd.Series(polyobject_pd.index).apply(lambda x: f'{polyline_object.metadata.name}_{x+1:04d}')
    
    #make gdf of points (1 point per row)
    gdf = geopandas.GeoDataFrame(data=df, geometry=geopandas.points_from_xy(polyobject_pd['x'],polyobject_pd['y']), crs=crs)
    
    return gdf


def PolyFile_to_geodataframe_points(polyfile_object:hcdfm.PolyFile, crs:str=None, add_pointnames:bool=True):
    """
    
    
    Parameters
    ----------
    polyfile_object : hcdfm.PolyFile
        get this object with hcdfm.PolyFile(path_to_plifile).
    crs : str, optional
        DESCRIPTION. The default is None'.
    add_pointnames : bool, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    gdf : TYPE
        gdf of all the pli files with columns: location, geometry (Points). Crs = 4326 (decimal degrees).

    """
    
    gdf_list = []
    for polyobj in polyfile_object.objects:
        gdf_one = pointlike_to_geodataframe_points(polyobj,crs=crs, add_pointnames=add_pointnames)
        gdf_list.append(gdf_one)
    gdf = pd.concat(gdf_list, ignore_index=True)
    return gdf


def da_from_gdf_points(gdf_points):
    
    ncbnd_construct = get_ncbnd_construct()
    dimn_point = ncbnd_construct['dimn_point']
    varn_pointx = ncbnd_construct['varn_pointx']
    varn_pointy = ncbnd_construct['varn_pointy']
    varn_pointname = ncbnd_construct['varn_pointname']
    attrs_pointx = ncbnd_construct['attrs_pointx']
    attrs_pointy = ncbnd_construct['attrs_pointy']
    
    da_plipoints = xr.Dataset()
    da_plipoints[varn_pointx] = xr.DataArray(gdf_points.geometry.x.tolist(), dims=dimn_point).assign_attrs(attrs_pointx)
    da_plipoints[varn_pointy] = xr.DataArray(gdf_points.geometry.y.tolist(), dims=dimn_point).assign_attrs(attrs_pointy)
    da_plipoints[varn_pointname] = xr.DataArray(gdf_points[varn_pointname].tolist(), dims=dimn_point)
    da_plipoints = da_plipoints.set_coords([varn_pointx,varn_pointy,varn_pointname])
    
    return da_plipoints


def PolyFile_to_geodataframe_linestrings(polyfile_object, crs=None):
    """
    empty docstring
    """
    plilines_list = []
    plinames_list = []
    pliz_list = []
    plidata_list = [] #TODO: make more generic with polyobject_pd.columns or earlier
    for iPO, polyline_object in enumerate(polyfile_object.objects):
        polyobject_pd = pd.DataFrame([dict(p) for p in polyline_object.points]) #TODO: getting only x/y might be faster, but maybe we also need the other columns like z/n/data?
        polygon_geom = LineString(zip(polyobject_pd['x'],polyobject_pd['y']))

        #make gdf of points (1 point per row)
        plilines_list.append(polygon_geom)
        plinames_list.append(polyline_object.metadata.name)
        pliz_list.append(polyobject_pd['z'].tolist())
        plidata_list.append(polyobject_pd['data'].tolist())

    gdf_polyfile = geopandas.GeoDataFrame({'name':plinames_list, 'z':pliz_list, 'data':plidata_list, 'geometry':plilines_list}, crs=crs)
    return gdf_polyfile


def gdf_linestrings_to_points(gdf_linestrings):

    crs = gdf_linestrings.crs
    gdf_points_list = []
    for ir, gdf_row in gdf_linestrings.iterrows():
        cols_list = list(gdf_linestrings.columns)
        cols_list.remove('geometry')

        gdf_one = geopandas.GeoDataFrame(data=None, geometry=geopandas.points_from_xy(*gdf_row.geometry.xy), crs=crs)
        for colname in cols_list:
            gdf_one[colname] = gdf_row[colname]
        gdf_points_list.append(gdf_one)

    gdf_points = pd.concat(gdf_points_list)
    gdf_points = gdf_points.reset_index(drop=True)
    return gdf_points


def parse_xy_to_datetime(pointlike_pd):
    datatimevals_pdstr = (pointlike_pd['x'].astype(int).apply(lambda x:f'{x:08d}') +
                          pointlike_pd['y'].astype(int).apply(lambda x:f'{x:06d}'))
    pointlike_pd.index = pd.to_datetime(datatimevals_pdstr)
    pointlike_pd_timeidx = pointlike_pd.drop(['x','y'],axis=1)
    
    return pointlike_pd_timeidx


def tekalobject_to_DataFrame(polyobject):
    #conversion to dataframe
    polyobject_pd = pointlike_to_DataFrame(polyobject)
    
    #conversion of xy to datetime
    polyobject_pd = parse_xy_to_datetime(polyobject_pd) #TODO: not suported the other way round, necessary to add?
    polyobject_pd.columns = polyobject.description.content.split('\n')[2:] #to fix legend labels
    polyobject_pd.index.name = polyobject.metadata.name
    return polyobject_pd


def maybe_convert_fews_to_dfmt(ds):
    ncbnd_construct = get_ncbnd_construct()
    dimn_point = ncbnd_construct['dimn_point']
    varn_pointname = ncbnd_construct['varn_pointname']
    
    # convert station data_vars to coords to avoid dfmt issues
    for var_to_coord in [varn_pointname,'station_names']:
        if var_to_coord in ds.data_vars:
            ds = ds.set_coords(var_to_coord)
    # assign timeseries_id cf_role to let FM read the station names
    ds[varn_pointname] = ds[varn_pointname].assign_attrs({'cf_role': 'timeseries_id'})
    
    # rename data_vars to long_name (e.g. renames FEWS so to salinitybnd)
    for datavar in ds.data_vars:
        if datavar in ['ux','uy']: #TODO: keeping these is consistent with hardcoded behaviour in dfm_tools elsewhere, but not desireable
            continue
        if hasattr(ds[datavar],'long_name'):
            longname = ds[datavar].attrs['long_name']
            if longname.endswith('bnd'):
                ds = ds.rename_vars({datavar:longname})
    
    # transpose dims #TODO: the order impacts the model results: https://issuetracker.deltares.nl/browse/UNST-7402
    # dfmt (arbitrary) dimension ordering is node/time/z
    # required to reorder to FEWS time/node/z order for comparable results
    # also time/z/node will result in unexpected results
    if "time" in ds.dims: # check if time dimension is present (astronomic does not have time)
        ds = ds.transpose("time", dimn_point, ...)
    
    # convert station names to string format (keep attrs and encoding)
    # also needed to properly export, since we cannot encode it at dtype S1 properly otherwise
    if not ds[varn_pointname].dtype.str.startswith('<'):
        with xr.set_options(keep_attrs=True):
            ds[varn_pointname] = ds[varn_pointname].load().str.decode('utf-8',errors='ignore').str.strip() #.load() is essential to convert not only first letter of string.

    # add relevant encoding if not present
    ds[varn_pointname].encoding.update({'dtype': 'S1', 'char_dim_name': 'char_leng_id'})
    
    # assign global attributes #TODO: add more
    ds = ds.assign_attrs({'Conventions': 'CF-1.6',
                          'institution': 'Deltares'})
    return ds
