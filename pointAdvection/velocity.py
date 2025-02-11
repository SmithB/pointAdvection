import pointCollection as pc
import geopandas as gpd
import os
import numpy as np
import shapely

class velocity(pc.grid.data):
    def __init__(self):
        self.t=None
        super(velocity, self).__init__()
        self.ndim=2

    # Purpose: Read velocity data from a file of known format
    def from_file(self, filename, format='NSIDC-0725', **kwargs):

        if format in ['NSIDC-0725','NSIDC-0727', 'NSIDC-0731', 'NSIDC-0766']:
            field_dict={'U':'vx','V':'vy','eU':'ex','eV':'ey'}
            self=self.from_h5(filename, field_mapping=field_dict, **kwargs)
        elif format=='NSIDC-0720' or format=='NSIDC-0484':
            # Antarctic ice velocity from UCI

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
