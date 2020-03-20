"""
Streamline plotting for 2D vector fields.
"""
from __future__ import division
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.collections as mcollections
import matplotlib.patches as patches
import bisect

__all__ = ['streamplot']


def streamplot(axes, x, y, u, v, density=1, linewidth=None, color=None,
               cmap=None, norm=None, arrowsize=1, arrowstyle='-|>',
               minlength=0.1, transform=None):
    """Draws streamlines of a vector flow.
    *x*, *y* : 1d arrays
        defines the grid.
    *u*, *v* : 2d arrays
        x and y-velocities. Number of rows should match length of y, and
        the number of columns should match x.
    *density* : float or 2-tuple
        Controls the closeness of streamlines. When `density = 1`, the domain
        is divided into a 25x25 grid---*density* linearly scales this grid.
        Each cell in the grid can have, at most, one traversing streamline.
        For different densities in each direction, use [density_x, density_y].
    *linewidth* : numeric or 2d array
        vary linewidth when given a 2d array with the same shape as velocities.
    *color* : matplotlib color code, or 2d array
        Streamline color. When given an array with the same shape as
        velocities, *color* values are converted to colors using *cmap*.
    *cmap* : :class:`~matplotlib.colors.Colormap`
        Colormap used to plot streamlines and arrows. Only necessary when using
        an array input for *color*.
    *norm* : :class:`~matplotlib.colors.Normalize`
        Normalize object used to scale luminance data to 0, 1. If None, stretch
        (min, max) to (0, 1). Only necessary when *color* is an array.
    *arrowsize* : float
        Factor scale arrow size.
    *arrowstyle* : str
        Arrow style specification.
        See :class:`~matplotlib.patches.FancyArrowPatch`.
    *minlength* : float
        Minimum length of streamline in axes coordinates.
    Returns:
        *stream_container* : StreamplotSet
            Container object with attributes
                - lines: `matplotlib.collections.LineCollection` of streamlines
                - arrows: collection of `matplotlib.patches.FancyArrowPatch`
                  objects representing arrows half-way along stream
                  lines.
            This container will probably change in the future to allow changes
            to the colormap, alpha, etc. for both lines and arrows, but these
            changes should be backward compatible.
    """
    grid = Grid(x, y)
    # Handle decreasing x and y by changing sign. The sign is changed
    # back after the integration routines (not totally happy with
    # this).
    x_increasing = grid.x[1] > grid.x[0]
    if not x_increasing:
        grid = Grid(-grid.x, grid.y)
        u = -u
    y_increasing = grid.y[1] > grid.y[0]
    if not y_increasing:
        grid = Grid(grid.x, -grid.y)
        v = -v

    mask = StreamMask(density)
    dmap = DomainMap(grid, mask)

    # default to data coordinates
    if transform is None:
        transform = axes.transData

    if color is None:
        color = axes._get_lines.color_cycle.next()

    if linewidth is None:
        linewidth = matplotlib.rcParams['lines.linewidth']

    line_kw = {}
    arrow_kw = dict(arrowstyle=arrowstyle, mutation_scale=10 * arrowsize)

    use_multicolor_lines = isinstance(color, np.ndarray)
    if use_multicolor_lines:
        assert color.shape == grid.shape
        line_colors = []
        if np.any(np.isnan(color)):
            color = np.ma.array(color, mask=np.isnan(color))
    else:
        line_kw['color'] = color
        arrow_kw['color'] = color

    if isinstance(linewidth, np.ndarray):
        assert linewidth.shape == grid.shape
        line_kw['linewidth'] = []
    else:
        line_kw['linewidth'] = linewidth
        arrow_kw['linewidth'] = linewidth

    ## Sanity checks.
    assert u.shape == grid.shape
    assert v.shape == grid.shape

    if np.any(np.isnan(u)):
        u = np.ma.array(u, mask=np.isnan(u))
    if np.any(np.isnan(v)):
        v = np.ma.array(v, mask=np.isnan(v))

    integrate = get_integrator(grid.x, grid.y, u, v, dmap, minlength)

    trajectories = []
    for xm, ym in _gen_starting_points(mask.shape):
        if mask[ym, xm] == 0:
            xg, yg = dmap.mask2data(xm, ym)
            t = integrate(xg, yg)
            if t is not None:
                trajectories.append(t)

    if use_multicolor_lines:
        if norm is None:
            norm = mcolors.normalize(color.min(), color.max())
        if cmap is None:
            cmap = cm.get_cmap(matplotlib.rcParams['image.cmap'])
        else:
            cmap = cm.get_cmap(cmap)

    streamlines = []
    arrows = []
    for t in trajectories:
        tx = np.array(t[0])
        ty = np.array(t[1])

        ## undoes the sign change put in place for the handling of
        ## decreasing x and y arrays.
        if not x_increasing:
            tx = -tx
        if not y_increasing:
            ty = -ty

        points = np.transpose([tx, ty]).reshape(-1, 1, 2)
        streamlines.extend(np.hstack([points[:-1], points[1:]]))

        # Add arrows half way along each trajectory.
        s = np.cumsum(np.sqrt(np.diff(tx) ** 2 + np.diff(ty) ** 2))
        n = np.searchsorted(s, s[-1] / 2.)
        arrow_tail = (tx[n], ty[n])
        arrow_head = (np.mean(tx[n:n + 2]), np.mean(ty[n:n + 2]))

        if isinstance(linewidth, np.ndarray):
            line_widths = interparray(grid, linewidth, tx, ty)[:-1]
            line_kw['linewidth'].extend(line_widths)
            arrow_kw['linewidth'] = line_widths[n]

        if use_multicolor_lines:
            color_values = interparray(grid, color, tx, ty)[:-1]
            line_colors.extend(color_values)
            arrow_kw['color'] = cmap(norm(color_values[n]))

        p = patches.FancyArrowPatch(arrow_tail,
                                    arrow_head,
                                    transform=transform,
                                    **arrow_kw)
        axes.add_patch(p)
        arrows.append(p)

    lc = mcollections.LineCollection(streamlines,
                                     transform=transform,
                                     **line_kw)
    if use_multicolor_lines:
        lc.set_array(np.asarray(line_colors))
        lc.set_cmap(cmap)
        lc.set_norm(norm)
    axes.add_collection(lc)

    axes.update_datalim(((x.min(), y.min()), (x.max(), y.max())))
    axes.autoscale_view(tight=True)

    ac = matplotlib.collections.PatchCollection(arrows)
    stream_container = StreamplotSet(lc, ac)
    return stream_container


class StreamplotSet(object):

    def __init__(self, lines, arrows, **kwargs):
        self.lines = lines
        self.arrows = arrows


# Coordinate definitions
#========================

class DomainMap(object):
    """Map representing different coordinate systems.
    Coordinate definitions:
    * axes-coordinates goes from 0 to 1 in the domain.
    * data-coordinates are specified by the input x-y coordinates.
    * mask-coordinates goes from 0 to N and 0 to M for an N x M mask,
      where N and M are user-specified to control the density of streamlines.
    This class also has methods for adding trajectories to the StreamMask.
    Before adding a trajectory, run `start_trajectory` to keep track of regions
    crossed by a given trajectory. Later, if you decide the trajectory is bad
    (e.g. if the trajectory is very short) just call `undo_trajectory`.
    """

    def __init__(self, grid, mask):
        self.grid = grid
        self.mask = mask
        ## Constants for conversion between grid- and mask-coordinates
        self.x_data2mask = float(mask.nx - 1) / grid.width
        self.y_data2mask = float(mask.ny - 1) / grid.height

        self.x_mask2data = 1. / self.x_data2mask
        self.y_mask2data = 1. / self.y_data2mask

    def data2mask(self, xi, yi):
        """Return nearest space in mask-coords from given data-coords."""
        return int((xi - self.grid.x_origin) * self.x_data2mask + 0.5), \
            int((yi - self.grid.y_origin)  * self.y_data2mask + 0.5)

    def mask2data(self, xm, ym):
        return self.grid.x_origin + xm * self.x_mask2data, \
            self.grid.y_origin + ym * self.y_mask2data

    def start_trajectory(self, xg, yg):
        xm, ym = self.data2mask(xg, yg)
        self.mask._start_trajectory(xm, ym)

    def reset_start_point(self, xg, yg):
        xm, ym = self.data2mask(xg, yg)
        self.mask._current_xy = (xm, ym)

    def update_trajectory(self, xg, yg):
        if not self.grid.within_grid(xg, yg):
            raise InvalidIndexError
        xm, ym = self.data2mask(xg, yg)
        self.mask._update_trajectory(xm, ym)

    def undo_trajectory(self):
        self.mask._undo_trajectory()


class Grid(object):
    """Grid of data."""
    def __init__(self, x, y):

        if len(x.shape) == 2:
            x_row = x[0]
            assert np.allclose(x_row, x)
            x = x_row
        else:
            assert len(x.shape) == 1

        if len(y.shape) == 2:
            y_col = y[:, 0]
            assert np.allclose(y_col, y.T)
            y = y_col
        else:
            assert len(y.shape) == 1

        self.nx = len(x)
        self.ny = len(y)

        self.x = x
        self.y = y

        self.x_origin = x[0]
        self.y_origin = y[0]
        self.width = x[-1] - x[0]
        self.height = y[-1] - y[0]

    @property
    def shape(self):
        return self.ny, self.nx

    def within_grid(self, xi, yi):
        """Return True if point is a valid index of grid."""
        return self.x_origin <= xi < self.x_origin + self.width and \
            self.y_origin <= yi < self.y_origin + self.height

class StreamMask(object):
    """Mask to keep track of discrete regions crossed by streamlines.
    The resolution of this grid determines the approximate spacing between
    trajectories. Streamlines are only allowed to pass through zeroed cells:
    When a streamline enters a cell, that cell is set to 1, and no new
    streamlines are allowed to enter.
    """

    def __init__(self, density):
        if np.isscalar(density):
            assert density > 0
            self.nx = self.ny = int(30 * density)
        else:
            assert len(density) == 2
            self.nx = int(25 * density[0])
            self.ny = int(25 * density[1])
        self._mask = np.zeros((self.ny, self.nx))
        self.shape = self._mask.shape

        self._current_xy = None

    def __getitem__(self, *args):
        return self._mask.__getitem__(*args)

    def _start_trajectory(self, xm, ym):
        """Start recording streamline trajectory"""
        self._traj = []
        self._update_trajectory(xm, ym)

    def _undo_trajectory(self):
        """Remove current trajectory from mask"""
        for t in self._traj:
            self._mask.__setitem__(t, 0)

    def _update_trajectory(self, xm, ym):
        """Update current trajectory position in mask.
        If the new position has already been filled, raise `InvalidIndexError`.
        """
        if self._current_xy != (xm, ym):
            if self[ym, xm] == 0:
                self._traj.append((ym, xm))
                self._mask[ym, xm] = 1
                self._current_xy = (xm, ym)
            else:
                raise InvalidIndexError


class InvalidIndexError(Exception):
    pass


class TerminateTrajectory(Exception):
    pass


# Integrator definitions
#========================

## This integrator now operates in *real space*.

def index_frac(x, x0):
    index = bisect.bisect(x, x0) - 1
    if index < 0: raise IndexError
    if index > len(x)-2: raise IndexError
    frac = (x0 - x[index]) / (x[index+1] - x[index])
    return index, frac


def get_integrator(x, y, u, v, dmap, minlength):

    # speed (path length) will be in axes-coordinates
    u_ax = u / dmap.grid.width
    v_ax = v / dmap.grid.height
    speed = np.ma.sqrt(u_ax ** 2 + v_ax ** 2)

    def forward_time(xi, yi):
        i_x, dx = index_frac(x, xi)
        i_y, dy = index_frac(y, yi)

        ds_dt = interpgrid(speed, i_x, dx, i_y, dy)
        if ds_dt == 0:
            raise TerminateTrajectory()
        dt_ds = 1. / ds_dt
        ui = interpgrid(u, i_x, dx, i_y, dy)
        vi = interpgrid(v, i_x, dx, i_y, dy)
        return ui * dt_ds, vi * dt_ds

    def backward_time(xi, yi):
        dxi, dyi = forward_time(xi, yi)
        return -dxi, -dyi

    def integrate(x0, y0):
        """Return x, y grid-coordinates of trajectory based on starting point.
        Integrate both forward and backward in time from starting point in
        grid coordinates.
        Integration is terminated when a trajectory reaches a domain boundary
        or when it crosses into an already occupied cell in the StreamMask. The
        resulting trajectory is None if it is shorter than `minlength`.
        """

        dmap.start_trajectory(x0, y0)
        sf, xf_traj, yf_traj = _integrate_rk12(x0, y0, dmap, forward_time)
        dmap.reset_start_point(x0, y0)
        sb, xb_traj, yb_traj = _integrate_rk12(x0, y0, dmap, backward_time)
        # combine forward and backward trajectories
        stotal = sf + sb
        x_traj = xb_traj[::-1] + xf_traj[1:]
        y_traj = yb_traj[::-1] + yf_traj[1:]

        if stotal > minlength:
            return x_traj, y_traj
        else:  # reject short trajectories
            dmap.undo_trajectory()
            return None

    return integrate


def _integrate_rk12(x0, y0, dmap, f):
    """2nd-order Runge-Kutta algorithm with adaptive step size.
    This method is also referred to as the improved Euler's method, or Heun's
    method. This method is favored over higher-order methods because:
    1. To get decent looking trajectories and to sample every mask cell
       on the trajectory we need a small timestep, so a lower order
       solver doesn't hurt us unless the data is *very* high resolution.
       In fact, for cases where the user inputs
       data smaller or of similar grid size to the mask grid, the higher
       order corrections are negligible because of the very fast linear
       interpolation used in `interpgrid`.
    2. For high resolution input data (i.e. beyond the mask
       resolution), we must reduce the timestep. Therefore, an adaptive
       timestep is more suited to the problem as this would be very hard
       to judge automatically otherwise.
    This integrator is about 1.5 - 2x as fast as both the RK4 and RK45
    solvers in most setups on my machine. I would recommend removing the
    other two to keep things simple.
    """
    ## This error is below that needed to match the RK4 integrator. It
    ## is set for visual reasons -- too low and corners start
    ## appearing ugly and jagged. Can be tuned.
    maxerror = 0.003

    ## This limit is important (for all integrators) to avoid the
    ## trajectory skipping some mask cells. We could relax this
    ## condition if we use the code which is commented out below to
    ## increment the location gradually. However, due to the efficient
    ## nature of the interpolation, this doesn't boost speed by much
    ## for quite a bit of complexity.
    maxds = min(1. / dmap.mask.nx, 1. / dmap.mask.ny, 0.1)

    ds = maxds
    stotal = 0
    xi = x0
    yi = y0
    xf_traj = []
    yf_traj = []

    while dmap.grid.within_grid(xi, yi):
        xf_traj.append(xi)
        yf_traj.append(yi)
        try:
            k1x, k1y = f(xi, yi)
            k2x, k2y = f(xi + ds * k1x,
                         yi + ds * k1y)
        except IndexError:
            # Out of the domain on one of the intermediate integration steps.
            # Take an Euler step to the boundary to improve neatness.
            ds, xf_traj, yf_traj = _euler_step(xf_traj, yf_traj, dmap, f)
            stotal += ds
            break
        except TerminateTrajectory:
            break

        dx1 = ds * k1x
        dy1 = ds * k1y
        dx2 = ds * 0.5 * (k1x + k2x)
        dy2 = ds * 0.5 * (k1y + k2y)

        # Error is normalized to the axes coordinates
        error = np.sqrt(((dx2 - dx1) / dmap.grid.width) ** 2
                        + ((dy2 - dy1) / dmap.grid.height) ** 2)

        # Only save step if within error tolerance
        if error < maxerror:
            xi += dx2
            yi += dy2
            try:
                dmap.update_trajectory(xi, yi)
            except InvalidIndexError:
                break
            if (stotal + ds) > 2:
                break
            stotal += ds

        # recalculate stepsize based on step error
        if error == 0:
            ds = maxds
        else:
            ds = min(maxds, 0.85 * ds * (maxerror / error) ** 0.5)

    return stotal, xf_traj, yf_traj


def _euler_step(xf_traj, yf_traj, dmap, f):
    """Simple Euler integration step that extends streamline to boundary."""
    ny, nx = dmap.grid.shape
    xi = xf_traj[-1]
    yi = yf_traj[-1]
    cx, cy = f(xi, yi)
    if cx == 0:
        dsx = np.inf
    elif cx < 0:
        dsx = (xi - dmap.grid.x_origin) / -cx
    else:
        dsx = (dmap.grid.x[-1] - xi) / cx
    if cy == 0:
        dsy = np.inf
    elif cy < 0:
        dsy = (yi - dmap.grid.y_origin) / -cy
    else:
        dsy = (dmap.grid.y[-1] - yi) / cy
    ds = min(dsx, dsy)
    xf_traj.append(xi + cx * ds)
    yf_traj.append(yi + cy * ds)
    return ds, xf_traj, yf_traj


# Utility functions
#========================
def interpgrid(a, i_x, dx, i_y, dy):
    """Fast 2D, linear interpolation on an integer grid.
    i_x, i_y are integer indices of array a corresponding to the cell of
           the array in which the point lies.
    dx, dy are the fractional indices corresponding to the location
           within the cell. These are typically between 0 and 1 (not
           required.)
    """

    a00 = a[i_y, i_x]
    a01 = a[i_y, i_x+1]
    a10 = a[i_y+1, i_x]
    a11 = a[i_y+1, i_x+1]
    a0 = a00 * (1 - dx) + a01 * dx
    a1 = a10 * (1 - dx) + a11 * dx
    ai = a0 * (1 - dy) + a1 * dy

    if not isinstance(i_x, np.ndarray):
        if np.ma.is_masked(ai):
            raise TerminateTrajectory

    return ai

def interparray(grid, a, x, y):
    """A simple 2d interpolation routine at points x and y (numpy
    arrays). Returns a numpy array of a at points x[i], y[i]."""

    i_x = np.clip((x[:,None] < grid.x).argmax(-1), 0, len(grid.x)-2)
    i_y = np.clip((y[:,None] < grid.y).argmax(-1), 0, len(grid.y)-2)
    i_x[x >= grid.x[-1]] = len(grid.x)-2
    i_y[y >= grid.y[-1]] = len(grid.y)-2

    dx = (x - grid.x[i_x]) / (grid.x[i_x+1] - grid.x[i_x])
    dy = (y - grid.y[i_y]) / (grid.y[i_y+1] - grid.y[i_y])
    return interpgrid(a, i_x, dx, i_y, dy)

def _gen_starting_points(shape):
    """Yield starting points for streamlines.
    Trying points on the boundary first gives higher quality streamlines.
    This algorithm starts with a point on the mask corner and spirals inward.
    This algorithm is inefficient, but fast compared to rest of streamplot.
    """
    ny, nx = shape
    xfirst = 0
    yfirst = 1
    xlast = nx - 1
    ylast = ny - 1
    x, y = 0, 0
    i = 0
    direction = 'right'
    for i in range(nx * ny):

        yield x, y

        if direction == 'right':
            x += 1
            if x >= xlast:
                xlast -= 1
                direction = 'up'
        elif direction == 'up':
            y += 1
            if y >= ylast:
                ylast -= 1
                direction = 'left'
        elif direction == 'left':
            x -= 1
            if x <= xfirst:
                xfirst += 1
                direction = 'down'
        elif direction == 'down':
            y -= 1
            if y <= yfirst:
                yfirst += 1
                direction = 'right'
