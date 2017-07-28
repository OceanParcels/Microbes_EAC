from parcels import (FieldSet, ParticleSet, JITParticle, AdvectionRK4, random,
                     ErrorCode, ParticleFile, Variable)
from datetime import timedelta as delta
import datetime
import numpy as np
import math


ddir = '/Volumes/data02/HYCOMdata/GLBa0.08_expt90_surf/hycom_GLBu0.08_912_'


def set_hycom_fieldset(t=None):
    t0 = datetime.date(2015, 6, 13)
    if t is None:
        files = [ddir + (t0 - delta(days=i)).strftime('%Y%m%d') + "00_t000.nc" for i in range(3, -1, -1)]
    else:
        files = [ddir + (t0 - delta(days=t)).strftime('%Y%m%d') + "00_t000.nc"]
    print files
    filenames = {'U': files, 'V': files, 'temp': files}
    variables = {'U': 'u', 'V': 'v', 'temp': 'temperature'}
    dimensions = {'lat': 'Latitude', 'lon': 'Longitude', 'time': 'MT', 'depth': 'Depth'}
    indices = {'lon': range(800, 1800), 'lat': range(700, 1600)}
    fset = FieldSet.from_netcdf(filenames, variables, dimensions, indices)
    fset.temp.data[fset.temp.data == 0] = np.nan
    return fset


def SampleTemp(particle, fieldset, time, dt):
    particle.temp = fieldset.temp[time, particle.lon, particle.lat, particle.depth]


def Age(particle, fieldset, time, dt):
    particle.age = particle.age + math.fabs(dt)
    if particle.age > fieldset.maxage:
        particle.delete()


def BrownianDiffusion(particle, fieldset, time, dt):
    kh_zonal = fieldset.Kh / math.pow(1000. * 1.852 * 60. * math.cos(particle.lat * M_PI / 180), 2)
    kh_meridional = fieldset.Kh / math.pow(1000.0 * 1.852 * 60.0, 2)

    r = 1/3.
    particle.lat += random.uniform(-1., 1.)*math.sqrt(2*math.fabs(dt)*kh_meridional/r)
    particle.lon += random.uniform(-1., 1.)*math.sqrt(2*math.fabs(dt)*kh_zonal/r)


def DeleteParticle(particle, fieldset, time, dt):
    particle.delete()


def run_microbes_eac(outfile):
    fieldset = set_hycom_fieldset()
    fieldset.add_constant('maxage', 90.*86400)
    fieldset.Kh = 100.  # diffusion constant

    i0 = len(fieldset.U.time)

    nperloc = 100
    sitelat38 = np.tile([-32.788], [nperloc])
    sitelon38 = np.tile([153.785], [nperloc])
    sitelat40 = np.tile([-30.621], [nperloc])
    sitelon40 = np.tile([153.371], [nperloc])

    class MicrobeParticle(JITParticle):
        temp = Variable('temp', dtype=np.float32, initial=fieldset.temp)
        age = Variable('age', dtype=np.float32, initial=0.)

    pset = {}
    pfile = {}
    pset[40] = ParticleSet(fieldset=fieldset, pclass=MicrobeParticle, lon=sitelon40,
                           lat=sitelat40, time=fieldset.U.time[-1])
    pfile[40] = ParticleFile(outfile+'40', pset[40])
    pfile[40].write(pset[40], pset[40][0].time)

    kernels = pset[40].Kernel(AdvectionRK4) + BrownianDiffusion + SampleTemp + Age

    for s in range(i0, 90, 1):
        if s is i0+1:
            pset[38] = ParticleSet(fieldset=fieldset, pclass=MicrobeParticle, lon=sitelon38,
                                   lat=sitelat38, time=fieldset.U.time[-1])
            pfile[38] = ParticleFile(outfile + '38', pset[38])
            pfile[38].write(pset[38], pset[38][0].time)
        for p in pset:
            pset[p].execute(kernels, starttime=pset[p][0].time, runtime=delta(days=1),
                            dt=delta(minutes=-5),
                            recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
            pfile[p].write(pset[p], pset[p][0].time)
        fieldset.advancetime(set_hycom_fieldset(s))


outfile = "microbes_eac_hycom_particles"
run_microbes_eac(outfile)
