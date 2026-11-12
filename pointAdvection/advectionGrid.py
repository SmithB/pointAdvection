#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
advectionGrid.py
Grids of parcel positions for moving points to and from a Lagrangian epoch

An advectionGrid has one or both pairs of these fields on a common (x, y, time) grid:
    x0, y0: the position at the epoch (time_0) of the parcel that is at (x, y) at time t
    x1, y1: the position at time t of the parcel that is at (x, y) at the epoch
Its time coordinate and time_0 are absolute, both in its units.

Files store x0, y0 in an 'xy0' group and x1, y1 in an 'xy1' group; a file
needs at least one of the two.  advectionGrid.to_h5 and to_nc write each
group with a 'time_0' attribute and absolute times, both in the units given
by the 'units' attribute of the time dataset (years if it is missing).

Older interpolator (xy01_grids) files, written by
scripts/make_vel_interpolator_for_xy.py, instead have a 'lagrangian_epoch'
group attribute (decimal years) and times relative to that epoch.
advectionGrid reads these and converts their times to absolute decimal years.

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    pointCollection: Utilities for organizing and manipulating point data
        https://github.com/SmithB/pointCollection
    timescale: Python tools for time and astronomical calculations
        https://github.com/pyTMD/timescale
"""
import warnings
import numpy as np
import pointAdvection.utilities

# attempt imports
h5py = pointAdvection.utilities.import_dependency('h5py')
netCDF4 = pointAdvection.utilities.import_dependency('netCDF4')
pc = pointAdvection.utilities.import_dependency('pointCollection')
timescale = pointAdvection.utilities.import_dependency('timescale')

try:
    obj = pc.grid.data
except AttributeError:
    obj = object

# the pair of position fields stored in each group of an interpolator file
_group_fields = {'xy0': ['x0', 'y0'], 'xy1': ['x1', 'y1']}
# the group that holds each position field
_field_group = {field: group for group, fields in _group_fields.items() for field in fields}
# descriptions written with the position fields
_descriptions = {
    'x0': 'x coordinate of the location at time_0 of a point observed at x, y at time t',
    'y0': 'y coordinate of the location at time_0 of a point observed at x, y at time t',
    'x1': 'x coordinate of the location at time t of a point that was at x, y at time_0',
    'y1': 'y coordinate of the location at time t of a point that was at x, y at time_0'}


def _unit_scale(from_units, to_units):
    """
    Factor that converts a time interval from from_units to to_units
    """
    try:
        return timescale.time._to_sec[from_units]/timescale.time._to_sec[to_units]
    except KeyError as exc:
        raise ValueError(f'unknown time units: {exc}') from exc


def _read_file_attributes(h5_file):
    """
    Find the position fields in a file, and read their epoch and time units

    Returns
    -------
    groups : list of str
        the groups ('xy0' and/or 'xy1') that contain both of their position fields
    time_0 : float
        the epoch: the 'time_0' attribute of the groups, in units, or for an
        older file the 'lagrangian_epoch' attribute, in decimal years
    units : str
        units of the times in the file (years if not specified)
    relative : bool
        True for an older file, whose times are relative to its lagrangian_epoch
    extent : list
        [xmin, xmax, ymin, ymax] of the grid
    """
    groups = []
    attrs = []
    extent = None
    with h5py.File(h5_file, 'r') as h5f:
        for group, fields in _group_fields.items():
            if group not in h5f or not all(field in h5f[group] for field in fields):
                continue
            g = h5f[group]
            if extent is None:
                x, y = np.array(g['x']), np.array(g['y'])
                extent = [np.min(x), np.max(x), np.min(y), np.max(y)]
            if 'time_0' in g.attrs:
                time_0, relative = float(np.squeeze(g.attrs['time_0'])), False
            elif 'lagrangian_epoch' in g.attrs:
                time_0, relative = float(np.squeeze(g.attrs['lagrangian_epoch'])), True
            else:
                raise ValueError(f'{h5_file}:{group} has no time_0 attribute')
            # pointCollection writes the time coordinate as 't' or 'time'
            units = 'years'
            for t_name in ('time', 't'):
                if t_name in g and 'units' in g[t_name].attrs:
                    units = g[t_name].attrs['units']
                    if isinstance(units, bytes):
                        units = units.decode('utf-8')
                    break
            groups.append(group)
            attrs.append((time_0, str(units), relative))
    if not groups:
        raise ValueError(f'{h5_file} contains neither x0 and y0 (in an xy0 group) '
                         'nor x1 and y1 (in an xy1 group)')
    if len(set(attrs)) > 1:
        raise ValueError(f'{h5_file}: the xy0 and xy1 groups have different epochs, '
                         'time units, or time formats')
    return (groups, *attrs[0], extent)


def _overlaps(extent, bounds):
    """
    Check whether a grid extent [xmin, xmax, ymin, ymax] overlaps bounds [[xmin, xmax], [ymin, ymax]]
    """
    if bounds is None:
        return True
    return (extent[0] <= bounds[0][1]) and (extent[1] >= bounds[0][0]) and \
        (extent[2] <= bounds[1][1]) and (extent[3] >= bounds[1][0])


class advectionGrid(obj):
    """
    Grid of parcel positions (x0, y0 and/or x1, y1), with absolute times

    Attributes
    ----------
    time_0: float or NoneType
        Lagrangian epoch to and from which parcels are advected
    units: str, default 'years'
        units of time_0 and of the time coordinate ('years' means decimal years)
    """
    def __init__(self, time_0=None, units='years', **kwargs):
        super().__init__(**kwargs)
        self.time_0 = time_0
        self.units = units

    def __copy__(self, fields=None):
        """Return a copy of the grid that keeps its class, time_0, and units."""
        return self.__from_data__(super().__copy__(fields=fields))

    def copy_meta(self):
        """Return an empty grid matching the current one."""
        return self.__from_data__(super().copy_meta())

    def __from_data__(self, temp):
        """
        Make an advectionGrid with this grid's time_0 and units from
        the attributes of a pointCollection.grid.data object
        """
        out = advectionGrid(time_0=self.time_0, units=self.units)
        out.__dict__.update(temp.__dict__)
        return out

    @property
    def groups(self):
        """
        The groups ('xy0' and/or 'xy1') whose position fields this grid has
        """
        return [group for group, fields in _group_fields.items()
                if all(field in self.fields for field in fields)]

    @staticmethod
    def read_time_0(files):
        """
        Read the epoch of advectionGrid or interpolator (xy01_grids) files
        without reading their grids

        Parameters
        ----------
        files: str or list
            hdf5 files with xy0 and/or xy1 groups

        Returns
        -------
        time_0: float
            the epoch of the files, in units
        units: str
            units of time_0 and of the times that from_h5 and from_files
            read ('years' means decimal years: older files, with times
            relative to their epoch, are read with absolute decimal years)
        """
        if isinstance(files, str):
            files = [files]
        if len(files) == 0:
            raise ValueError('advectionGrid.read_time_0: no files')
        epochs = {}
        for file in files:
            _, time_0, units, relative, _ = _read_file_attributes(file)
            epochs[file] = (time_0, 'years' if relative else units)
        if len(set(epochs.values())) > 1:
            raise ValueError('interpolator files have different epochs or time units: ' +
                             ', '.join(f'{file}: time_0={t0} {units}' for file, (t0, units) in epochs.items()))
        return epochs[files[0]]

    def __merge__(self, grids, time_0, units, relative, copy=False):
        """
        Set self from grids holding x0, y0 and/or x1, y1

        Parameters
        ----------
        grids: dict
            pointCollection grids, keyed by group ('xy0', 'xy1')
        time_0: float
            the epoch: in units if relative is False, otherwise in decimal years
        units: str
            units of the grids' times
        relative: bool
            if True, the grids' times are relative to time_0 and are
            converted to absolute decimal years
        """
        if not grids:
            raise ValueError('an advectionGrid needs x0 and y0, or x1 and y1')
        for group, grid in grids.items():
            missing = [field for field in _group_fields[group] if getattr(grid, field, None) is None]
            if missing:
                raise ValueError(f'the {group} grid has no {" or ".join(missing)} field')
        grid_list = list(grids.values())
        for other in grid_list[1:]:
            for coord in ('x', 'y', 'time'):
                c0, c1 = getattr(grid_list[0], coord), getattr(other, coord)
                if (c0 is None) != (c1 is None) or (c0 is not None and not np.array_equal(c0, c1)):
                    raise ValueError(f'the xy0 and xy1 grids have different {coord} coordinates')
        # take the coordinates and metadata from the first grid, then add the position fields
        self.__dict__.update(pc.grid.data.copy_meta(grid_list[0]).__dict__)
        for group, grid in grids.items():
            for field in _group_fields[group]:
                z = getattr(grid, field)
                self.assign({field: z.copy() if copy else z})
        self.__update_size_and_shape__()
        self.time_0 = time_0
        if not relative:
            self.units = units
            return self
        # convert the relative times to absolute decimal years
        self.units = 'years'
        if self.time is not None:
            self.time = time_0 + np.asarray(self.time, dtype=float)*_unit_scale(units, self.units)
        return self

    def __relative_t_range__(self, t_range, time_0, units):
        """
        Convert an absolute time range (decimal years) to times relative to time_0 in units
        """
        if t_range is None:
            return None
        return [(t - time_0)*_unit_scale('years', units) for t in t_range]

    def from_h5(self, h5_file, t_range=None, **kwargs):
        """
        Read an advectionGrid or interpolator (xy01_grids) file

        Parameters
        ----------
        h5_file: str
            hdf5 file with an xy0 group (x0, y0), an xy1 group (x1, y1), or both
        t_range: list or NoneType, default None
            absolute time range to read, in the file's time units
            (decimal years for an older file with relative times)
        kwargs: dict
            other keywords for pointCollection.grid.data.from_h5 (e.g. bounds)
        """
        groups, time_0, units, relative, _ = _read_file_attributes(h5_file)
        if relative:
            t_range = self.__relative_t_range__(t_range, time_0, units)
        grids = {group: pc.grid.data().from_h5(h5_file, group=group, fields=_group_fields[group],
                                               t_range=t_range, **kwargs)
                 for group in groups}
        return self.__merge__(grids, time_0, units, relative)

    def from_files(self, files, bounds=None, **kwargs):
        """
        Mosaic a list of advectionGrid or interpolator (xy01_grids) files

        The mosaic includes the position fields (x0, y0 and/or x1, y1) that
        every file overlapping the bounds has.

        Parameters
        ----------
        files: list
            hdf5 files with xy0 and/or xy1 groups, all with the same epoch,
            time units, and time format (absolute, or relative for older files)
        bounds: list or NoneType, default None
            boundaries to read: ``[[xmin, xmax], [ymin, ymax]]``
        kwargs: dict
            other keywords for pointCollection.grid.mosaic.from_list
        """
        if isinstance(files, str):
            files = [files]
        info = [_read_file_attributes(file) for file in files]
        if len(set(file_info[1:4] for file_info in info)) > 1:
            details = ', '.join(f'{file}: epoch={t0}, units={units}, '
                                f'{"relative" if relative else "absolute"} times'
                for file, (_, t0, units, relative, _) in zip(files, info))
            raise ValueError('interpolator files have different epochs, time units, '
                             f'or time formats: {details}')
        _, time_0, units, relative, _ = info[0]
        # only the files that overlap the bounds go into the mosaic
        used = [(file, file_info[0]) for file, file_info in zip(files, info)
                if _overlaps(file_info[-1], bounds)]
        if not used:
            raise ValueError(f'none of the interpolator files overlap the bounds {bounds}: ' + ', '.join(files))
        files = [file for file, _ in used]
        groups = []
        for group in _group_fields:
            lacking = [file for file, file_groups in used if group not in file_groups]
            if not lacking:
                groups.append(group)
            elif len(lacking) < len(used):
                warnings.warn(f'advectionGrid: not reading {" and ".join(_group_fields[group])}, '
                              f'which are missing from {", ".join(lacking)}')
        if not groups:
            raise ValueError('the interpolator files do not all contain x0 and y0, '
                             'and do not all contain x1 and y1: ' + ', '.join(files))
        grids = {}
        for group in groups:
            # mosaic.from_list removes files outside the bounds from its input list: pass a copy
            temp = pc.grid.mosaic().from_list(list(files), group=group, fields=_group_fields[group],
                bounds=bounds, **kwargs)
            if isinstance(temp, str):
                # mosaic.from_list returns a message when it cannot read the inputs
                raise ValueError(temp)
            grids[group] = temp
        return self.__merge__(grids, time_0, units, relative)

    def from_grid(self, xy0=None, xy1=None, time_0=None, units='years'):
        """
        Make an advectionGrid from grids with times relative to time_0,
        such as the outputs of advection.xy0_interpolator and xy1_interpolator

        Parameters
        ----------
        xy0, xy1: pointCollection.grid.data or NoneType
            grids with fields x0, y0 and x1, y1.  At least one is needed, and
            if both are given they must have the same x, y, and time coordinates
        time_0: float
            Lagrangian epoch, in decimal years
        units: str, default 'years'
            units of the relative times in xy0 and xy1
        """
        if time_0 is None:
            raise ValueError('time_0 must be specified')
        grids = {group: grid for group, grid in (('xy0', xy0), ('xy1', xy1)) if grid is not None}
        return self.__merge__(grids, time_0, units, relative=True, copy=True)

    def interp(self, *args, **kwargs):
        """
        Interpolate a field (see pointCollection.grid.data.interp)

        Times are absolute, in self.units.
        """
        field = kwargs.get('field', 'z')
        if field in _field_group and field not in self.fields:
            raise ValueError(f'this advectionGrid has no {field} field: '
                             f'its input had no {_field_group[field]} group')
        return super().interp(*args, **kwargs)

    def __write_groups__(self, writer, out_file, replace, **kwargs):
        """
        Write each group of position fields with writer (the pointCollection
        to_h5 or to_nc method)

        Returns
        -------
        groups: list
            the groups that were written
        """
        if self.time_0 is None:
            raise ValueError('time_0 must be set to write an advectionGrid')
        groups = self.groups
        if not groups:
            raise ValueError('an advectionGrid needs x0 and y0, or x1 and y1')
        for group in groups:
            writer(out_file, group=group, fields=_group_fields[group], replace=replace, **kwargs)
            # only the first group replaces the file
            replace = False
        return groups

    def to_h5(self, out_file, replace=False, **kwargs):
        """
        Write the grid to an hdf5 file

        x0, y0 are written to the 'xy0' group and x1, y1 to the 'xy1' group
        (for whichever pairs the grid has), with absolute times.  Each group
        has a 'time_0' attribute, the time dataset has a 'units' attribute
        (the units of time_0 and of the times), and the position datasets
        have 'description' attributes.

        Parameters
        ----------
        out_file: str
            hdf5 file
        replace: bool, default False
            overwrite the file (otherwise the groups are added to it)
        kwargs: dict
            other keywords for pointCollection.grid.data.to_h5
        """
        groups = self.__write_groups__(super().to_h5, out_file, replace, **kwargs)
        with h5py.File(out_file, 'r+') as h5f:
            for group in groups:
                g = h5f[group]
                g.attrs['time_0'] = self.time_0
                for t_name in ('time', 't'):
                    if t_name in g:
                        g[t_name].attrs['units'] = self.units
                for field in _group_fields[group]:
                    g[field].attrs['description'] = _descriptions[field]
        return self

    def to_nc(self, out_file, replace=False, **kwargs):
        """
        Write the grid to a netCDF4 file

        The layout is the same as for to_h5: x0, y0 are written to an 'xy0'
        group and x1, y1 to an 'xy1' group (for whichever pairs the grid has),
        with absolute times.  Each group has a 'time_0' attribute, the time
        variable has a 'units' attribute (the units of time_0 and of the
        times), and the position variables have 'description' attributes.

        Parameters
        ----------
        out_file: str
            netCDF4 file
        replace: bool, default False
            overwrite the file (otherwise the groups are added to it)
        kwargs: dict
            other keywords for pointCollection.grid.data.to_nc
        """
        kwargs['t_units'] = self.units
        groups = self.__write_groups__(super().to_nc, out_file, replace, **kwargs)
        with netCDF4.Dataset(out_file, 'a') as ncf:
            for group in groups:
                g = ncf.groups[group]
                g.setncattr('time_0', self.time_0)
                for field in _group_fields[group]:
                    g.variables[field].setncattr('description', _descriptions[field])
        return self
