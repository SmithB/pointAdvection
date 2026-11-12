#!/usr/bin/env python3
"""
advectionGrid_for_tile.py
Make an advectionGrid (x0, y0, x1, y1 on an x, y, time grid) for a tile
from a set of velocity maps

This does the same calculation as make_vel_interpolator_for_xy.py, and takes
the same arguments, but writes the result as an advectionGrid: x0, y0 in an
'xy0' group and x1, y1 in an 'xy1' group, with absolute times (decimal years)
and the Lagrangian epoch in each group's 'time_0' attribute.  The output is
read with pointAdvection.advectionGrid().from_h5 or from_files.

The values are the same as make_vel_interpolator_for_xy.py's, but they are
calculated faster: only the nodes inside the tile are tracked (unless
--no_crop is set), --processes splits the tracking between processes, and
the velocity interpolation uses numba if it is installed.

Example:
    advectionGrid_for_tile.py --xy0 2400000 -400000 @default_args/Shack_vel_interpolator_args.txt --processes 8
"""
import os
import sys
import json
import argparse
import multiprocessing
import numpy as np
import pointCollection as pc
import pointAdvection


def parse_range(value):
    """
    Parse a range given as one comma-separated string or as two strings
    """
    if value is None:
        return None
    if ',' in value[0]:
        return [*map(float, value[0].split(','))]
    return [*map(float, value)]


def parse_input_args(argv):
    """
    Transform input arguments into a namespace

    Parameters
    ----------
    argv : iterable
        Input arguments.  Keywords should have two hyphens at the start;
        @file reads arguments from a file

    Returns
    -------
    argparse.Namespace
        input arguments
    """
    parser = argparse.ArgumentParser(
        description='Make an advectionGrid (displacement grids) for a tile '
                    'based on a set of velocity measurements.',
        fromfile_prefix_chars='@')
    parser.add_argument('--xy0', type=float, nargs=2, required=True, help='tile center')
    parser.add_argument('--width', '-w', type=float, required=True, help='tile width')
    parser.add_argument('--w_smooth', type=float,
                        help='fill gaps in the mean velocity field (a "mean" json entry) '
                             'using this smoothing length')
    parser.add_argument('--t_range', type=str, nargs='+', required=True,
                        help='starting and ending times: two numbers comma or space separated')
    parser.add_argument('--velocity_t_range', type=str, nargs='+',
                        help='range of velocity measurements to use: two numbers comma or space separated')
    parser.add_argument('--velocity_source_file', type=str, required=True,
                        help='json file listing times, filenames, and formats, '
                             'or a velocity file with a time dimension (see --velocity_format)')
    parser.add_argument('--velocity_format', type=str, default='NSIDC-0731',
                        help='format of velocity_source_file if it is not a json file: NSIDC-0731, '
                             'or h5 for a velocity stack with x, y, time, U, V, eU, and eV datasets')
    parser.add_argument('--res', type=float, default=500,
                        help='not used: the output grid has the resolution of the velocity maps '
                             '(kept so that make_vel_interpolator_for_xy.py argument files work)')
    parser.add_argument('--lagrangian_epoch', type=float, required=True,
                        help='time (time_0) to/from which points will be advected')
    parser.add_argument('--t_step', type=float, default=0.25, help='time resolution of the output grids')
    parser.add_argument('--advection_time_step', type=float, default=1/48,
                        help='time step for parcel tracking, in years')
    parser.add_argument('--speed_est', type=float, default=1.e3,
                        help='speed estimate used to pad the area of velocity data read')
    parser.add_argument('--save_base', type=str, default='./xy0_interpolators')
    parser.add_argument('--EPSG', type=int, default=3031, help='SRS EPSG for the velocity maps')
    parser.add_argument('--no_crop', action='store_true',
                        help='unless this is set, the output will be cropped to the specified width')
    parser.add_argument('--quiet', action='store_true', help='no verbose output')
    parser.add_argument('--interpolator_save_file', type=str,
                        help='output file (.h5, or .nc for netCDF4)')
    parser.add_argument('--processes', type=int, default=1,
                        help='number of processes for the parcel tracking (the results do not depend on this)')

    args, _ = parser.parse_known_args(argv)
    args.VERBOSE = not args.quiet
    args.t_range = parse_range(args.t_range)
    args.velocity_t_range = parse_range(args.velocity_t_range)
    return args


def make_adv_from_json(args, xg, yg, bds_pad):
    """
    Make an advection object from velocity files listed in a json file

    Entries with a numeric 'time' are velocity maps; an entry with
    time 'mean' is a mean velocity map used to fill gaps.
    """
    with open(args.velocity_source_file, 'r') as fh:
        vel_files = json.load(fh)

    entries = [vi for vi in vel_files if not isinstance(vi['time'], str)]
    if args.velocity_t_range is not None:
        entries = [vi for vi in entries
                   if args.velocity_t_range[0] <= vi['time'] <= args.velocity_t_range[1]]
    entries = sorted(entries, key=lambda vi: vi['time'])

    # make a list of advection objects, one for each velocity file
    advs = []
    for vi in entries:
        temp = pointAdvection.velocity().from_file(vi['file'], format=vi['format'], bounds=bds_pad)
        if 'mask' not in temp.fields:
            temp.assign(mask=np.isfinite(temp.U).astype(float))
        if temp.time is None:
            temp.time = vi['time']
        adv_i = pointAdvection.advection()
        adv_i.x = xg
        adv_i.y = yg
        adv_i.t = vi['time']
        adv_i.velocity = temp
        advs += [adv_i]

    adv = pointAdvection.advection().from_list(advs)
    # if there's a mean velocity file, read it
    for vi in vel_files:
        if vi['time'] == 'mean':
            adv.vel_mean = pointAdvection.velocity().from_file(vi['file'], format=vi['format'], bounds=bds_pad)
            if args.w_smooth is not None:
                N_smooth = np.ceil(args.w_smooth/(adv.vel_mean.x[1]-adv.vel_mean.x[0]))
                adv.vel_mean.fill_smoothed(fields=['U', 'V', 'eU', 'eV'], w_smooth=N_smooth)
            adv.vel_mean = adv.vel_mean.interp_to(adv.velocity.x, adv.velocity.y)
    return adv


def make_adv_from_stack(args, xg, yg, bds_pad):
    """
    Make an advection object from a velocity file with a time dimension:
    an NSIDC-0731 (GrIMP) file, or (with --velocity_format=h5) a velocity
    stack such as those written by pointCollection.grid.data.to_h5
    """
    adv = pointAdvection.advection()
    adv.x = xg
    adv.y = yg
    adv.velocity = pointAdvection.velocity().from_file(args.velocity_source_file,
                                                       format=args.velocity_format, bounds=bds_pad)
    adv.t0 = args.lagrangian_epoch
    adv.velocity.spacing = (adv.velocity.x[1]-adv.velocity.x[0], adv.velocity.y[1]-adv.velocity.y[0])
    adv.velocity.type = 'grid'
    return adv


def default_save_file(args):
    """
    Output filename used by make_vel_interpolator_for_xy.py
    """
    return f'{args.save_base}/E{int(args.xy0[0]/1000)}_N{int(args.xy0[1]/1000)}_' + \
        f'{args.t_range[0]}_{args.lagrangian_epoch}_{args.t_range[1]}.h5'


def prepare_advection(args):
    """
    Read the velocity maps for a tile, fill their gaps, and set them up for
    parcel tracking: times relative to the epoch, and velocities scaled to
    projected-coordinate displacements

    Returns
    -------
    pointAdvection.advection
    """
    bds = [jj + np.array([-1, 1])*args.width*0.5 for jj in args.xy0]
    pad = args.speed_est*np.diff(args.t_range)*np.array([-1, 1])*2
    bds_pad = [jj+pad for jj in bds]
    xg, yg = np.meshgrid(*[np.arange(bdsi[0], bdsi[1], args.res) for bdsi in bds])

    if args.velocity_source_file.endswith('.json'):
        adv = make_adv_from_json(args, xg, yg, bds_pad)
    else:
        adv = make_adv_from_stack(args, xg, yg, bds_pad)

    # if there was no 'mean' velocity file, calculate the mean velocity
    if adv.vel_mean is None:
        adv.vel_mean = adv.calc_vel_mean(w_smooth=5)

    # fill velocity gaps
    adv.fill_velocity_gaps(annual=False, vel_mean=adv.vel_mean)

    # the advection calculations use times relative to the epoch
    t_vel = adv.velocity.time if adv.velocity.time is not None else adv.velocity.t
    adv.velocity.time = np.array(t_vel, dtype=float) - args.lagrangian_epoch

    # account for the ratio between velocity and projected-coordinate displacements
    # (ps_scale_for_lat is only correct for the standard polar stereographic grids)
    if args.EPSG in (3031, 3413):
        adv.velocity.get_latlon(srs_epsg=args.EPSG)
        scale = pc.ps_scale_for_lat(adv.velocity.latitude)
        if adv.velocity.U.ndim == 3:
            scale = scale[:, :, None]
        adv.velocity.U /= scale
        adv.velocity.V /= scale
    elif args.VERBOSE:
        print(f'EPSG={args.EPSG} is not 3031 or 3413: '
              'velocities will not be corrected for projection distortion')
    adv.t0 = 0
    return adv


def crop_bounds(args):
    """
    x and y ranges of the tile (the output is cropped to these)
    """
    return [args.xy0[0]+args.width*np.array([-0.5, 0.5]),
            args.xy0[1]+args.width*np.array([-0.5, 0.5])]


def interpolator_bounds(velocity, bounds=None):
    """
    Bounds to pass to xy0_interpolator and xy1_interpolator

    The interpolators put their grid nodes at b0-dx, b0, b0+dx, ... up to
    b1, for bounds [b0, b1] in x and y, with dx the x spacing of the
    velocity grid.  The default is the bounds of the velocity grid.  Every
    node is tracked independently, so to make the nodes inside a smaller
    area, it is enough to pass bounds that are on the same lattice: this
    gives the same nodes, with the same values, inside that area.

    Parameters
    ----------
    velocity: pointCollection.grid.data
        velocity grid
    bounds: list or NoneType
        x and y ranges that the nodes should cover, or None for the whole
        velocity grid
    """
    full = velocity.bounds()
    if bounds is None:
        return full
    dx = velocity.x[1] - velocity.x[0]
    out = []
    for (b0, b1), (lo, hi) in zip(full, bounds):
        k0 = np.floor((lo - b0)/dx)
        k1 = np.ceil((hi - b0)/dx)
        out.append([max(b0 + k0*dx, b0), min(b0 + k1*dx, b1)])
    return out


# advection object used by worker processes (inherited when they are forked)
_worker_adv = None


def _interpolators_for_bounds(bounds, kwargs, blocksize=200):
    """
    Make the xy1 and xy0 grids for one set of bounds
    """
    adv = _worker_adv
    # pass separate copies of the bounds (older versions of the interpolators appended t_range to them).
    # nodes are tracked independently, so the block size affects only the speed
    # (blocks of 200x200 nodes were faster than make_vel_interpolator_for_xy.py's 100 and 50)
    xy1 = adv.xy1_interpolator(bounds=[list(b) for b in bounds], blocksize=blocksize, **kwargs)
    xy0 = adv.xy0_interpolator(bounds=[list(b) for b in bounds], blocksize=blocksize, **kwargs)
    return xy1, xy0


def run_interpolators(adv, bounds, processes=1, chunk=None, **kwargs):
    """
    Make the xy1 and xy0 grids, optionally in parallel

    With more than one process, the grid is split into chunks of
    chunk x chunk nodes that are tracked in separate processes and then
    joined: each node is tracked independently, so the result is the same
    as with one process.

    Parameters
    ----------
    adv: pointAdvection.advection
        advection object with the velocity grid
    bounds: list
        bounds for the interpolators (see interpolator_bounds)
    processes: int
        number of processes
    chunk: int or NoneType
        size of the chunks of nodes given to the processes.  The default
        gives about two chunks per process, with chunks of 50 to 100 nodes
        (smaller chunks are slowed by per-call overhead)
    kwargs: dict
        t_range, t_step, and advection_time_step for the interpolators
    """
    global _worker_adv
    _worker_adv = adv
    if processes <= 1:
        return _interpolators_for_bounds(bounds, kwargs)
    # the nodes that the interpolators make for these bounds
    dx = adv.velocity.x[1] - adv.velocity.x[0]
    nodes = [np.arange(bb[0]-dx, bb[1]+dx, dx) for bb in bounds]
    if chunk is None:
        chunk = int(np.clip(np.ceil(np.sqrt(nodes[0].size*nodes[1].size/(2*processes))), 50, 100))
    # chunks of node indices in x and y, and bounds that make the nodes in each chunk
    splits = [[ind[i:i+chunk] for i in range(0, ind.size, chunk)]
              for ind in [np.arange(n.size) for n in nodes]]
    chunks = [(cols, rows) for rows in splits[1] for cols in splits[0]]
    tasks = [([[nodes[0][cols[0]]+dx, nodes[0][cols[-1]]], [nodes[1][rows[0]]+dx, nodes[1][rows[-1]]]],
              kwargs, None) for cols, rows in chunks]
    # build the interpolation tables before forking so that the processes share them
    adv.__linear_grid__()
    with multiprocessing.get_context('fork').Pool(processes) as pool:
        results = pool.starmap(_interpolators_for_bounds, tasks, chunksize=1)
    joined = []
    for ii, fields in enumerate([['x1', 'y1'], ['x0', 'y0']]):
        t = results[0][ii].t
        out = {'x': nodes[0], 'y': nodes[1], 't': t}
        for field in fields:
            out[field] = np.full((nodes[1].size, nodes[0].size, t.size), np.nan)
        for (cols, rows), result in zip(chunks, results):
            grid = result[ii]
            # place each chunk's own nodes (rounding might add a node at a chunk's edge)
            keep_x = np.abs(grid.x[:, None] - nodes[0][cols][None, :]).min(axis=1) < dx/2
            keep_y = np.abs(grid.y[:, None] - nodes[1][rows][None, :]).min(axis=1) < dx/2
            for field in fields:
                out[field][rows[0]:rows[-1]+1, cols[0]:cols[-1]+1, :] = \
                    getattr(grid, field)[np.ix_(keep_y, keep_x)]
        joined.append(pc.grid.data().from_dict(out))
    return tuple(joined)


def advectionGrid_for_tile(argv):
    """
    Make and save the advectionGrid for a tile

    Returns
    -------
    pointAdvection.advectionGrid
    """
    args = parse_input_args(argv)
    if args.interpolator_save_file is not None:
        out_file = args.interpolator_save_file
    else:
        out_file = default_save_file(args)
    if args.VERBOSE:
        print(f"working on {out_file}")

    adv = prepare_advection(args)
    # when the output is cropped, track only the nodes that it keeps
    bounds = interpolator_bounds(adv.velocity, None if args.no_crop else crop_bounds(args))
    xy1, xy0 = run_interpolators(adv, bounds, processes=args.processes,
                                 t_range=np.array(args.t_range) - args.lagrangian_epoch,
                                 t_step=args.t_step, advection_time_step=args.advection_time_step)
    if not args.no_crop:
        for grid in [xy0, xy1]:
            grid.crop(*crop_bounds(args))

    out = pointAdvection.advectionGrid().from_grid(xy0=xy0, xy1=xy1, time_0=args.lagrangian_epoch,
                                                   units='years')
    if os.path.dirname(out_file):
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
    if out_file.endswith('.nc'):
        out.to_nc(out_file, replace=True)
    else:
        out.to_h5(out_file, replace=True)
    return out


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    advectionGrid_for_tile(argv)


if __name__ == '__main__':
    main()
