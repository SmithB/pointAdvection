import pointAdvection.utilities

# attempt imports
pc = pointAdvection.utilities.import_dependency('pointCollection')

try:
    obj = pc.grid.data
except AttributeError:
    obj = object

import geopandas as gpd
import os
import numpy as np
import shapely
import re
import datetime
class velocity(obj):
    def __init__(self):
        self.t=None
        super(velocity, self).__init__()
        self.ndim=2


    def _from_local_vrt(self, v_file, read_errors=True, **kwargs):


        t_re = re.compile('Vel-(\d\d\d\d)-(\d\d)-(\d\d).(\d\d\d\d)-(\d\d)-(\d\d).*')
        y0, m0, d0, y1, m1, d1 = map(int, t_re.search(v_file).groups())
        t0=datetime.datetime(y0, m0, d0)
        t1=datetime.datetime(y1, m1, d1)
        t = t0 + (t1-t0)/2

        self.from_geotif(v_file, meta_only=True, **kwargs)
        field_dict={'U':'vx','V':'vy','eU':'ex','eV':'ey'}
        if read_errors:
            field_dict['eU']='ex'
            field_dict['eV']='ey'
        for out_field, in_field in field_dict.items():
            temp=pc.grid.data().from_file(v_file.replace('v.vrt',in_field+'.vrt'), **kwargs)
            self.assign({out_field:temp.z})
        self.t = 2018 + (t-datetime.datetime(2018, 1, 1)).days/365.25


    def _from_NSIDC_GL(self, filename, **kwargs):
        field_dict={'U':'vx','V':'vy','eU':'ex','eV':'ey'}
        self.from_h5(filename, field_mapping=field_dict, **kwargs)

    def _from_NSIDC_AA(self, filename, **kwargs):
        field_dict={'U':'VX','V':'VY','eU':'ERRX','eV':'ERRY'}
        self=self.from_nc(filename, field_mapping=field_dict, **kwargs)
        self.EPSG=3031
        edits_file=filename.replace('.nc','_edits.geojson')
        if os.path.isfile(edits_file):
            self.assign(mask=np.ones_like(self.U))
            df=gpd.read_file(edits_file)
            for poly in df.geometry:
                self.rasterize_poly(poly, field='mask', burn_value=0, raster_epsg=self.EPSG, poly_epsg=self.EPSG)
            if np.any(self.mask==0):
                self.mask[self.mask==0]=np.nan
                for field in self.fields:
                    if field=='mask':
                        continue
                    if len(self.U.shape)==2:
                        setattr(self, field, getattr(self, field) * self.mask)
                    else:
                        for band in self.U.shape[3]:
                            temp=getattr(self, field)[:,:,band]
                            temp *= self.mask
                            getattr(self, field)[:,:,band]=temp
    # Purpose: Read velocity data from a file of known format
    def from_file(self, filename, format='NSIDC-0725', **kwargs):

        if format in ['NSIDC-0725','NSIDC-0727', 'NSIDC-0731', 'NSIDC-0766']:
            self._from_NSIDC_GL(filename, **kwargs)
        elif format=='NSIDC-0720' or format=='NSIDC-0484':
            # Antarctic ice velocity from UCI
            self._from_NSIDC_AA(filename,**kwargs)
        elif format=='local_vrt':
            self._from_local_vrt(filename, **kwargs)
        else:
            print(f"format {format} for file {filename} unknown, skipping")
        self.ndim=len(self.U.shape)
        return self

    def interp_to(self, x0, y0, t0=None):

        temp={'x':x0,
            'y':y0,
            'time':t0}
        for field in self.fields:
            temp[field]=self.interp(x0.ravel(), y0.ravel(), t0, field=field, gridded=True)
        return velocity().from_dict(temp)
