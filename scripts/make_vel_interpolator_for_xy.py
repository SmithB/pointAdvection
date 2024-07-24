import numpy as np
import re
import h5py
import os
import pointAdvection
import json
import datetime
import re
import argparse
import sys

# sample input arguments:


def parse_input_args(argv):
    '''
    parse_input_args: transform input argument string into a dataspace

    Parameters
    ----------
    args : iterable
        Input arguments.  Keywords should have two hyphens at the start

    Returns
    -------
    dataspace
        input arguments formatted as a namespace

    '''

    parser = argparse.ArgumentParser(description='Make a displacement grid for a location based on a set of velocity measurements.', \
                                     fromfile_prefix_chars='@')
    parser.add_argument('--xy0', type=float, nargs=2, help='interpolator center')
    parser.add_argument('--width','-w', type=float, help='width of interpolator')
    parser.add_argument('--t_range', nargs=2, type=float, help='starting and ending times')
    parser.add_argument('--velocity_t_range', nargs=2, type=float, help='range of velocity measurements to use')


    parser.add_argument('--velocity_source_file', type=str, help='json file listing times, filenames, and formats')
    parser.add_argument('--res', type=float, default=500, help='resolution of output map')
    parser.add_argument('--lagrangian_epoch', type=float,  help='time to/from which points will be advected')

    parser.add_argument('--t_step', type=float, default=0.25, help='time resolution of output maps')
    parser.add_argument('--speed_est', type=float, default=1.e3, help='speed estimate used to pad input range')
    parser.add_argument('--save_base', type=str, default='./xy0_interpolators')
    parser.add_argument('--interpolator_save_file', type=str)

    args, unk=parser.parse_known_args(argv)
    return args

def make_adv_from_json(args, xg, yg):
    with open(args.velocity_source_file,'r') as fh:
        vel_files=json.load(fh)
 
    times = [vi['time'] for vi in vel_files if not isinstance(vi['time'],str)]
    if args.velocity_time_range is not None:
        times= [tt for tt in times if (tt >= args.velocity_time_range[0]) and (tt <= args.velocity_time_range[1])]

    time_order = np.argsort(times)

    # make a list of advection objects, one for each velocity file
    advs=[]

    for ind in time_order:
        temp=pointAdvection.velocity().from_file(vel_files[ind]['file'], format=vel_files[ind]['format'],  bounds=bds_pad)
        if temp.time is None:
            if temp.t is not None:
                temp.time=temp.t
            else:
                temp.time=vel_files[ind]['time']
        adv_i=pointAdvection.advection()
        adv_i.x=xg
        adv_i.y=yg
        adv_i.t=times[ind]
        adv_i.velocity=temp
        advs += [adv_i]

    adv = pointAdvection.advection().from_list(advs)    

def make_adv_from_grimp(args):
    temp = pointAdvection.velocity().from_file(args.velocity_source_file, format='NSIDC-0731',  bounds=bds_pad)    

    
def make_vel_interpolator_for_xy(argv):
    args=parse_input_args(argv)

    if args.interpolator_save_file is not None:
        interpolator_save_file=args.interpolator_save_file
    else:
        interpolator_save_file=f'{args.save_base}_{args.t_range[0]}_{args.lagrangian_epoch}_{args.t_range[1]}.h5'
    print(interpolator_save_file)

    
    bds = [jj + np.array([-1, 1])*args.width*0.5 for jj in args.xy0]

    pad=args.speed_est*np.diff(args.t_range)*np.array([-1, 1])
    bds_pad = [jj+pad for jj in bds]
    xg, yg = np.meshgrid(*[np.arange(bdsi[0], bdsi[1], args.res) for bdsi in bds])
    
    if args.velocity_source_file.endswith('.json'):
        adv = make_adv_from_json(args, xg, yg)
    else:
        adv = make_adv_from_grimp(args, xg, yg)
        
    # check if there is a 'mean' velocity file
    for vi in vel_files:
        if vi['time']=='mean':
          adv.vel_mean=pointAdvection.velocity().from_file(vi['file'], format=vi['format'], bounds=bds_pad)\
            .interp_to(adv.velocity.x, adv.velocity.y)
    if adv.vel_mean is None:
        adv.vel_mean=adv.calc_vel_mean(w_smooth=5)

    # fill velocity gaps
    adv.fill_velocity_gaps(annual=False, vel_mean=adv.vel_mean)


    # make the advection interpolators

    t_vel = adv.velocity.time.copy()
    adv.velocity.time = t_vel- args.lagrangian_epoch
    blocksize=100
    bounds=adv.velocity.bounds()
    # create an xy0 interpolator to find final locations for each point.
    xy1 = adv.xy1_interpolator(bounds=bounds,
                                   t_range=(np.array(args.t_range)-args.lagrangian_epoch),
                                   t_step=args.t_step, advection_time_step=1/12, blocksize=blocksize)
    # create an xy0 interpolator to find final locations for each point
    xy0=adv.xy0_interpolator(bounds=bounds,
                               t_range=(np.array(args.t_range)-args.lagrangian_epoch),
                               t_step=args.t_step, advection_time_step=1/12, blocksize=50)
    if os.path.isfile(interpolator_save_file):
        os.remove(interpolator_save_file)
    xy1.to_h5(interpolator_save_file, group='xy1', replace=True)
    xy0.to_h5(interpolator_save_file, group='xy0', replace=False)
    with h5py.File(interpolator_save_file,'r+') as h5f:
        for group in ['xy0','xy1']:
            h5f[group].attrs['lagrangian_epoch']=args.lagrangian_epoch
            h5f[group]['t'].attrs['units']='years'

def main(args=None):
    if args is None:
        args=sys.argv
    make_vel_interpolator_for_xy(args)

if __name__=='__main__':
    main()
