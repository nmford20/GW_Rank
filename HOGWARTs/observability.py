import astropy
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astroplan import Observer, FixedTarget, AtNightConstraint, AltitudeConstraint, is_observable

def check_obs(site,ras,decs):

     #select subset of galaxies visible at telescope location in next 24 hours (nighttime only)

     telescope=Observer.at_site(site)

     # Get time now

     T0 = astropy.time.Time.now()

     #get time range of next 24 hours

     time_range = [T0,T0+24*u.hour]

     #get the lists of ra and decs as columns of array

     coords = np.column_stack((ras,decs))

     print("grabbed coordinate pairs")
     #convert these coord pairs to targets

     targets = [FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg))
        for ra, dec in coords]

     print("converted coords to targets") 

     #constrain to nightime (using astronomical twilight) in the next 24 hours
     #And constrain to objects at least 8 deg above horizon (*CFHT constraint, modify for telescope*)

     constraints = [AtNightConstraint.twilight_astronomical(),AltitudeConstraint(8*u.deg, 90*u.deg)]

     # Are targets ever observable in the chosen time range?
     #Note: this will overestimate the amount of galaxies you can actually observe,
     #since it takes time to observe each galaxy and it doesn't account for how long/when a galaxy is observable

     ever_observable = np.array(is_observable(constraints, telescope, targets, time_range=time_range))

     print("found which targets are observable")

     #get the row numbers of targets that are observable

     idxs = list(np.where(ever_observable == True)[0])

     print("grabbed relevant indices")

     return idxs