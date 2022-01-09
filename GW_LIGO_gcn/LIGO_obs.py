
import gcn
import lxml.etree
import subprocess
import healpy as hp
from ligo.skymap.io import read_sky_map
import astropy
from astropy.io import fits
import urllib.request
import pathlib
import numpy as np
import astropy.coordinates
import astropy.time
import astropy.units as u
from GW_GalRank import GW_GalRank
import createfiles

"""
Provide telescope conditions here:
Should be a list in the form telescope geodetic coordinates lat (decimal deg), lon (decimal deg), height (m), 
sun altaz, airmass (secant of zenith angle approximation)
"""
Conditions = [19.825252,-155.468876,4204,-12,2.5]

#Checking GW observability
#Note: LIGO/Virgo probability sky maps are always in equatorial coordinates

"""
This function takes the following inputs: 
skymap, header, telescope geodetic coordinates latitude (deg), longitude (deg), height (m), angle of the sun, airmass,
    T0 (default = None, uses current time. If not None, takes a string in one of the astropy Time formats.)
"""

def prob_observable(skymap, header, la, lo, h, sun, airmass, T0=None):
    """
    Determine the integrated probability contained in a gravitational-wave
    sky map that is observable from a particular ground-based site at a
    particular time.

    Bonus: make a plot of probability versus UTC time!
    """

    if T0 is not None:
        #read in the astropy-formatted string containing the desired time, ex.'2015-03-01 13:55:27'
        time = astropy.time.Time(T0)
    else:
        # Get time now
        time = astropy.time.Time.now()
    # Or at the time of the gravitational-wave event...
    # time = astropy.time.Time(header['MJD-OBS'], format='mjd')

    # Geodetic coordinates of observatory
    observatory = astropy.coordinates.EarthLocation(
        lat=la*u.deg, lon=lo*u.deg, height=h*u.m)

    # Alt/az reference frame at observatory, now
    frame = astropy.coordinates.AltAz(obstime=time, location=observatory)

    # Look up (celestial) spherical polar coordinates of HEALPix grid.
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    # Convert to RA, Dec.
    radecs = astropy.coordinates.SkyCoord(
        ra=phi*u.rad, dec=(0.5*np.pi - theta)*u.rad)

    # Transform grid to alt/az coordinates at observatory, now
    altaz = radecs.transform_to(frame)

    # Where is the sun, now?
    sun_altaz = astropy.coordinates.get_sun(time).transform_to(altaz)

    # How likely is it that the (true, unknown) location of the source
    # is within the area that is visible, now? Demand that sun is at
    # least 18 degrees below the horizon and that the airmass
    # (secant of zenith angle approximation) is at most 2.5.
    prob = skymap[(sun_altaz.alt <= sun*u.deg) & (altaz.secz <= airmass)].sum()

    # Done!
    return prob

#Convenience function for calculating the time between GW signal and current time
#Values to be used for Tobs parameter in the tiling code
#talert = signal emission time (in Astropy Time units)

def tdiff(talert):
    #turning talert into an astropy Time object in case it isn't already in that form
    talert = astropy.time.Time(talert, format = 'gps')
    tcurrent = astropy.time.Time(astropy.time.Time.now().gps, format = 'gps')
    T0 = (tcurrent - talert).value
    #add 24 hours to get a 24 hour window for observations (tiling code will select nightime only)
    T1 = (tcurrent - talert + 24*u.hour).value
    return T0, T1

# Handler Function to call every time a GCN is received.
# Run only for notices of type
# LVC_PRELIMINARY, LVC_INITIAL, LVC_UPDATE, or LVC_RETRACTION.

@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE,
    gcn.notice_types.LVC_RETRACTION)
def process_gcn(payload, root):
    # Respond only to 'test' events.
    # VERY IMPORTANT! Replace with the following code
    # to respond to only real 'observation' events.
    # if root.attrib['role'] != 'observation':
    #    return
    if root.attrib['role'] != 'observation':
        return

    # Read all of the VOEvent parameters from the "What" section.
    #This will store the values into a params dictionary
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}

    if params['AlertType'] == 'Retraction':
        print(params['GraceID'], 'was retracted')
        return

    # Respond only to 'CBC' events. Change 'CBC' to 'Burst'
    # to respond to only unmodeled burst events.
    if params['Group'] != 'CBC':
        return

    graceid=params['GraceID']
    prelim=params['AlertType']+params['Pkt_Ser_Num']

    for key, value in params.items():
        print(key, '=', value)

    #If the alert includes a probability skymap, lets analyze it
    if 'skymap_fits' in params:

        #download a copy of the skymap to the 'skymaps' directory 
        #(create directory if needed)
        
        url_main = params['skymap_fits']
        file = pathlib.PurePath(url_main)
        
        #get the current working directory for downloading files/saving outputs
        cwd = str(pathlib.Path.cwd())

        local_path = str(cwd+'/../skymaps/'+graceid)
        filename = str(local_path+"/"+file.name)

        #check if any parent directories need to be created
        pathlib.Path(filename).parents[0].mkdir(exist_ok=True, parents=True)

        #save the file
        urllib.request.urlretrieve(url_main, filename)
        
        # Read the FITS header.
        hdulist = fits.open(filename)
        distest = hdulist[1].header['DISTMEAN']
        diststd= hdulist[1].header['DISTSTD']
        utctime = hdulist[1].header['DATE-OBS']
        gpstime = astropy.time.Time(utctime, scale='utc').gps

        # Print some values from the FITS header.
        print('Distance =', distest, '+/-', diststd)

        #read the HEALPIX skymap
        skymap, meta=read_sky_map(filename, moc=False,distances=True)

        #extracting some key information
        prob=skymap[0]
        distmu=skymap[1]
        distsigma=skymap[2]
        distnorm=skymap[3]
        npix = len(prob)
        nside=hp.npix2nside(npix)

        #OPTIONAL (uncomment below): Do an observability estimate, using pre-defined telescope conditions
        #prob0 = prob_observable(skymap, header, Conditions[0], Conditions[1], Conditions[2], Conditions[3], Conditions[4])
        #print('Source has a {:d}% chance of being observable now'.format(int(round(100 * prob0))))
        #if prob0 > 0.5:
            #pass # FIXME: perform some action
            
        print("extracted skymap parameters, moving on to galaxy ranking...")

        #running the HOGWARTs galaxy ranking code

        #Set where you want output list to go
        outputdir = cwd+"/.."+"/outputs/GW_GalRank_output/"+graceid+"/"
        #check if all the subfolders exist - if not, make them
        pathlib.Path(outputdir).parent.mkdir(exist_ok=True, parents=True)
        info, levelsper, ncontour = GW_GalRank(skymap,prob,distest,diststd,distmu,distsigma,distnorm,nside,graceid,prelim,url_main,local_path,outputdir)
        createfiles.createasciifile(info,graceid,prelim,levelsper,ncontour,outputdir)
        print("saved ranked galaxies list as ascii file, moving on to tile ranking...")

        #running the gwemopt tile ranking code
        T0, T1 = tdiff(gpstime)
        #need to feed in the shell command to execute, w/command line arguments
        subprocess.run("python gwemopt_run --telescope CFHT --do3D --doTiles --doPlots --doSchedule --timeallocationType powerlaw --scheduleType greedy -o '../../outputs/GW_TileRank_output/"+graceid+"/powerlaw_greedy' --gpstime "+str(float(gpstime))+" --doSkymap --tilesType ranked --skymap '../../skymaps/"+graceid+"/"+str(file.name)+"' --filters g --exposuretimes 1000.0 --Tobs "+str(T0)+","+str(T1)+" --doSingleExposure",shell=True,cwd="../gwemopt/bin")
        print("*COMPLETE* waiting for next gcn alert...")

#sample skymap to test script, note that you need to download the xml file first
#uncomment these three lines if you already downloaded the skymap, otherwise comment them out
payload = open('../skymaps/S190814bv/S190814bv-5-Update.xml', 'rb').read()
root = lxml.etree.fromstring(payload)
process_gcn(payload, root)

# Listen for GCNs until the program is interrupted
# (killed or interrupted with control-C).
gcn.listen(handler=process_gcn)