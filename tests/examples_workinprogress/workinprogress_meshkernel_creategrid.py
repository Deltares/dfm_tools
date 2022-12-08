# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 20:54:39 2022

@author: veenstra

"""

import meshkernel

import matplotlib.pyplot as plt
plt.close()
import numpy as np



def create_mk_instance_with_a_uniform_curvilinear_grid(num_columns =3, num_rows =3):
    r"""A local function for creating an instance of meshkernel with a uniform curvilinear grid.
    """
    mk = meshkernel.MeshKernel()
    
    # Create an instance of MakeGridParameters and set the values
    make_grid_parameters = meshkernel.MakeGridParameters() #TODO: does not exist in release code
    make_grid_parameters.num_columns = num_columns
    make_grid_parameters.num_rows = num_rows
    make_grid_parameters.angle = 0.0
    make_grid_parameters.block_size = 0.0
    make_grid_parameters.origin_x = 0.0
    make_grid_parameters.origin_y = 0.0
    make_grid_parameters.block_size_x = 10.0
    make_grid_parameters.block_size_y = 10.0
    
    # A polygon must to be provided. If empty it will not be used.
    node_x = np.empty(0, dtype=np.double)
    node_y = np.empty(0, dtype=np.double)
    geometry_list = meshkernel.GeometryList(node_x, node_y)

    mk.curvilinear_make_uniform(make_grid_parameters, geometry_list)
    
    return mk

## Make a uniform rectilinear grid
"""
1.	Make a regular and rotated grid. This is possible in MeshKernelPy but I need to merge the related pull request in master. The same workflow of interactor is followed here: first generate a curvilinear grid than convert the curvilinear grid into unstructured grid. The steps are the following:
a.	curvilinear_make_uniform, see the following notebook: https://github.com/Deltares/MeshKernelPy/blob/AddCurvilinearGridSupport/docs/examples/04_curvilineargrid_basics.ipynb
b.	curvilinear_convert_to_mesh2d: https://github.com/Deltares/MeshKernelPy/blob/118cb4953c4e95d5b18ed283bb37f391134b2bb2/meshkernel/meshkernel.py#L1399 

    If you want to write the resulting unstructured grid to file, you can use UGridPy. See the following example notebook: https://github.com/Deltares/UGridPy/blob/main/docs/examples/05_working_with_meshkernel.ipynb

"""
mk = create_mk_instance_with_a_uniform_curvilinear_grid()
curvilinear_grid = mk.curvilineargrid_get()
fig, ax = plt.subplots()
curvilinear_grid.plot_edges(ax)

#If a polygon is provided it will be used in the generation of the curvilinear grid. The polygon must be closed
node_x = np.array([2.5,5.5,3.5,0.5,2.5], dtype=np.double)
node_y = np.array([0.5,3.0,5.0,2.5,0.5], dtype=np.double)
geometry_list = meshkernel.GeometryList(node_x, node_y)

make_grid_parameters = meshkernel.MakeGridParameters()
make_grid_parameters.num_columns = 10
make_grid_parameters.num_rows = 10
make_grid_parameters.angle = 0.0
make_grid_parameters.block_size = 0.0
make_grid_parameters.origin_x = 0.0
make_grid_parameters.origin_y = 0.0
make_grid_parameters.block_size_x = 0.2
make_grid_parameters.block_size_y = 0.2

mk.curvilinear_make_uniform(make_grid_parameters, geometry_list)
curvilinear_grid = mk.curvilineargrid_get()
fig, ax = plt.subplots()
curvilinear_grid.plot_edges(ax)

#convert to ugrid/mesh2d
mk_ugrid = mk.curvilinear_convert_to_mesh() #TODO: does not exist in release code: https://github.com/Deltares/MeshKernelPy/blob/118cb4953c4e95d5b18ed283bb37f391134b2bb2/meshkernel/meshkernel.py#L1399




# refine with samples. #TODO: this does not seem to be related to bathymetry/courant refinement
"""
3.	Mesh refinement in MeshKernelPy (where the samples parameter contains the sample coordinates and the sample depths used for mesh refinement):

https://github.com/Deltares/MeshKernelPy/blob/6a5110d75ea2a86a27f63da7c7196e008bbce3dc/meshkernel/meshkernel.py#L451

"""
mk.mesh2d_refine_based_on_samples() #'samples', 'relative_search_radius', 'minimum_num_samples', and 'mesh_refinement_params'





