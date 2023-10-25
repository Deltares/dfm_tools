from dfm_tools import modelbuilder as mb

netfile = r'p:\11208614-de-370a\01_models\BES\DFLOWFM\Models\bonaire\bonaire_net.nc'

interval = 0.05
resolution = 'h'
threshold_mindepth = 20

obs = mb.generate_coastline_obspoints(interval = interval,
                                        resolution = resolution,
                                        threshold_mindepth = threshold_mindepth,
                                        file_nc = netfile)


