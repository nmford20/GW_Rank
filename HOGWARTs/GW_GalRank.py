import sys
import lxml.etree
import gcn
import healpy as hp
import pathlib
import urllib
import numpy as np
import astropy_healpix
import astropy.utils.data
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
import pandas as pd
from ligo.skymap.postprocess import crossmatch
from ligo.skymap.io import read_sky_map
import ligo.skymap.plot
from probability import calculate_absolute_probability,unique_galaxies,extract_LIGO_probability,makelists,sortbyprob
from dataframe import create_dataframe
from contours import contour_plots,join_0_360,split_contours,integrated_probability,hpix_contours
from createfiles import createtxt,createjsonfile,createasciifile
from GLADE import Gcoordinates
from skymapio import readskymap
from astroplan import Observer, FixedTarget, AtNightConstraint, AltitudeConstraint, is_observable
from observability import check_obs

"""
    Compute ranked list of galaxy probabilities using a ligo probability skymap.

    Parameters:
    -----------
    skymap: ligo HEALPIX skymap, extracted from ligo.skymap.io.read_sky_map(),

    prob: list 
        list of LIGO/Virgo galaxy probabilities 

    distest: float
        mean estimated distance to GW event

    diststd: float
        standard deviation in mean estimated distance

    distmu: list 
        list of skymap mean distances

    distssigma: list 
        list of skymap distance standard deviations

    distnorm: list 
        list of skymap normalisation factors 

    nside: int
        resolution of skymap

    graceid: string
        Assigned identifier for GW event

    skymap_path: string (optional)
        directory path to skymap file. Default to local 'skymaps' folder

    output_path: string (optional)
        directory path to save output contour plot and final ranking list. Default to local 
        'GW_GalRank_output' subfolder named by event's graceid

    Return:
    -------
    Outputs a saved ascii file containing the final ranked galaxy list
    in descending probability order, w/coordinates
    
    """

def GalRank(skymap,prob,distest,diststd,distmu,distsigma,distnorm,nside,graceid,prelim,url_main,local_path,outputdir):

    #create integrated probability skymap
    csm=integrated_probability(prob)
    #calculate contours from csm
    contours=hpix_contours(csm,levels=[0.90],nest=False)

    #define 90% region
    levels=[0.90]
    levelsper=[90]

    #define distance limits
    distmax=distest+5*diststd
    distmin=distest-5*diststd

    #get coordinates of all GLADE+ galaxies (w/some initial selection cuts)
    coordinates,data=Gcoordinates(distmax,distmin)

    print("Read GLADE+ catalog and get all galaxy coordinates")

    #crossmatch GLADE with multiorder skymap
    if 'v1' in url_main:
        version='.v1'
    elif 'v0' in url_main:
        version='.v0'
    elif 'v2' in url_main:
        version='.v2'
    else:
        version=''
    if ',0' in url_main:
        after=',0'
    elif ',1' in url_main:
        after=',1'
    elif ',2' in url_main:
        after=',2'
    elif ',3' in url_main:
        after=',3'
    else:
        after=''

    if 'bayestar' in url_main:
        url='https://gracedb.ligo.org/api/superevents/'+graceid+'/files/bayestar.multiorder.fits'+after
    else:
        url='https://gracedb.ligo.org/api/superevents/'+graceid+'/files/LALInference'+version+'.multiorder.fits'+after
    #note: encountered one skymap where the file had .offline. in the name, not sure why. It throws an error though
    file = pathlib.PurePath(url)
    urllib.request.urlretrieve(url, local_path+"/"+file.name)
    skymap_moc=read_sky_map(local_path+"/"+file.name, moc=True)
    result=crossmatch(skymap_moc,coordinates)

    #for each contour region (Eg. 99%)
    for d in range(0,len(contours)):

        #jsonlist=[]
        jsonlist2=[]
        tablenames=[]
        ra_incontourlist=[]
        contourlens=[]
        dec_incontourlist=[]
        finalprobslist=[]
        finalgalnamelist=[]
        dist_incontourlist=[]
        Bmag_incontourlist=[]
        ra_incontourlist1=[]
        dec_incontourlist1=[]
        probs_incontourlist1=[]
        probs_incontourlist=[]
        finalgalnamelist1=[]
        dist_incontourlist1=[]
        Bmag_incontourlist1=[]
        finalgalname=[]
        mudists_incontourlist1=[]
        distssigma_incontourlist1=[]
        distsnorm_incontourlist1=[]
        pdist_incontourlist1=[]
        Slum_incontourlist1=[]
        contourlist1=[]
        contourlist=[]
        contourss=[]
        ccc=[]
        
        #separate masked array into separate contours
        split_dec, split_ra = split_contours(contours, levels[d],d)

        #retrieve galaxies in 90 percent regions
        results=data[result.searched_prob<0.90]
        ra_incontour=results['RA'].values
        dec_incontour=results['Dec'].values
        dist_incontour=results['dist'].values
        Bmag_incontour=results['Bmag'].values
        name_incontour=results['HyperLEDA'].values
        
        # if the contour is split at 0/360 degrees, rejoin back together for plot
        split_ra2, split_dec2=join_0_360(split_ra, split_dec)
        
        #create a plot of contours and number them
        contour_plots(split_ra2, split_dec2,graceid, prelim, levelsper[d], outputdir)
        contourss=np.ones(len(ra_incontour))

        # extract probability parameters at galaxy positions
       
        probs, mudists, distssigma, distsnorm = extract_LIGO_probability(ra_incontour, dec_incontour, nside, distsigma, prob, distnorm, distmu)
        
        # remove duplicates
        indices, finalgalname, ra_incontourlist1,dec_incontourlist1,dist_incontourlist1,Bmag_incontourlist1,probs_incontourlist1,mudists_incontourlist1,distssigma_incontourlist1,distsnorm_incontourlist1, contourlist=unique_galaxies(contourlist, contourss,ra_incontourlist1, ra_incontour, dec_incontourlist1, dec_incontour, dist_incontourlist1, dist_incontour, probs_incontourlist1, name_incontour, finalgalname, probs,Bmag_incontourlist1, Bmag_incontour,mudists_incontourlist1, mudists,distssigma_incontourlist1, distssigma,distsnorm_incontourlist1, distsnorm)
        
        # Calculate probability score
       
        finalprobs,pdist,Slum = calculate_absolute_probability(dist_incontourlist1, Bmag_incontourlist1, mudists_incontourlist1, distssigma_incontourlist1,distsnorm_incontourlist1,probs_incontourlist1)
        
        finalprobss=[]
        for j in range(0,len(finalprobs[0])):
            finalprobss.append(finalprobs[0,j])
        # make lists for dataframes

        finalprobss,ra_incontourlist,dec_incontourlist,finalprobslist,probs_incontourlist, finalgalnamelist, dist_incontourlist,Bmag_incontourlist,contourlist1,pdist_incontourlist1, Slum_incontourlist1 = makelists(finalprobss,ra_incontourlist,ra_incontourlist1,dec_incontourlist,dec_incontourlist1,finalprobslist, probs_incontourlist1, probs_incontourlist, finalgalnamelist,finalgalname,dist_incontourlist,dist_incontourlist1,Bmag_incontourlist,Bmag_incontourlist1,contourlist1,contourlist,pdist, pdist_incontourlist1, Slum, Slum_incontourlist1)
        
        #sort by descending probability
                                                                                         
        finaldictsorted, cumsumprobs=sortbyprob(finalprobslist)
        
        print("sorted galaxies in descending probability order")

                
        #select subset of galaxies observable from telescope location in next 24 hours (nighttime only)

        idxs = check_obs('cfht',ra_incontourlist,dec_incontourlist)

        #create dataframe for jsons
        dataf=create_dataframe(finaldictsorted, ra_incontourlist, dec_incontourlist, probs_incontourlist, finalgalnamelist, dist_incontourlist, pdist_incontourlist1, Bmag_incontourlist, Slum_incontourlist1, contourlist, cumsumprobs)

        #selecting the indices we want
        dataf = dataf.iloc[idxs]

        #jsonlist.append(dataf[['Galaxy name', 'Galaxy probability score', 'RA (degrees)', 'Dec (degrees)','Location probability score','Distance (Mpc)', 'Distance probability score', 'B magnitude', 'B luminosity probability score','Cumulative Score']].to_json())
        jsonlist2.append(dataf[['Galaxy name', 'Galaxy probability score', 'RA (degrees)', 'Dec (degrees)','Location probability score','Distance (Mpc)', 'Distance probability score','B magnitude','B luminosity probability score','Cumulative Score']].to_csv())
        
        #createtxt(dataf,finalgalnamelist, finaldictsorted,graceid,prelim,levelsper,d,ccc)                                                                
        #createjsonfile(jsonlist,graceid,prelim,levelsper,d)
        return jsonlist2, levelsper, d
