
import os, sys
import copy
import numpy as np
import healpy as hp
import itertools

from scipy.stats import norm

import astropy.coordinates
from astropy.table import unique, vstack, Table
from astropy.time import Time, TimeDelta
import astropy.units as u


def read_transients(params, map_struct):

    transients = Table.read(params["transientsFile"], format="ascii",
                            names=("name","ra","dec"))

    nside = params["nside"]
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5*np.pi - theta)

    ipix = hp.ang2pix(nside, theta=transients["ra"], phi=transients["dec"],
                      lonlat=True)
    prob = map_struct["prob"][ipix]
    prob = prob / np.sum(prob)

    transients_struct = {}
    transients_struct["ra"] = transients["ra"]
    transients_struct["dec"] = transients["dec"]
    transients_struct["name"] = transients["name"]
    transients_struct["S"] = prob
    transients_struct["Sloc"] = prob
    transients_struct["Smass"] = prob

    return transients_struct

def combine_catalog_transients(params, catalog_struct, transients_struct):

    transient_weight = params["transients_to_catalog"]
    catalog_weight = 1 - transient_weight
    catalog_struct["ra"] = np.hstack([catalog_struct["ra"],
                                      transients_struct["ra"]])
    catalog_struct["dec"] = np.hstack([catalog_struct["dec"],
                                       transients_struct["dec"]])

    catalog_struct["Sloc"] = np.hstack([catalog_struct["Sloc"]*catalog_weight,
                                        transients_struct["Sloc"]*transient_weight])
    catalog_struct["S"] = np.hstack([catalog_struct["S"]*catalog_weight,
                                     transients_struct["S"]*transient_weight])
    catalog_struct["Smass"] = np.hstack([catalog_struct["Smass"]*catalog_weight,
                                         transients_struct["Smass"]*transient_weight])

    return catalog_struct

def read_ps1_transients(params, map_struct):

    nside = params["nside"]
    prob_data_sorted = np.sort(map_struct["prob"])[::-1]
    prob_data_indexes = np.argsort(map_struct["prob"])[::-1]
    prob_data_cumsum = np.cumsum(prob_data_sorted)

    lines = [line.rstrip('\n') for line in open(params["transientsFile"])]
    lines = lines[1:]
    lines = filter(None,lines)

    transients_struct = {}
    transients_struct["data"] = np.empty((0,8))
    transients_struct["name"] = []
    transients_struct["classification"] = []

    for line in lines:
        line = line.replace('"','')
        lineSplit = line.split(",")
        if lineSplit[0] == "NULL":
            continue

        name = lineSplit[0]
        ra = float(lineSplit[2])
        dec = float(lineSplit[3])
        classification = lineSplit[6]

        if not lineSplit[11] == "":
            mag = float(lineSplit[11])
        else:
            mag = -1

        if not lineSplit[13] == "":
            z = float(lineSplit[13])
        else:
            z = -1
        if not lineSplit[14] == "":
            mpc = float(lineSplit[14])
        else:
            mpc = -1

        if not lineSplit[15] == "":
            discovery_mjd = Time(lineSplit[15], format='isot').mjd
        else:
            discovery_mjd = -1

        if not lineSplit[16] == "":
            nondetection_mjd = Time(lineSplit[16], format='isot').mjd
        else:
            #nondetection_mjd = -1
            nondetection_mjd = discovery_mjd

        event_mjd = Time(params["gpstime"], format='gps', scale='utc').mjd
        if (event_mjd < nondetection_mjd) or (event_mjd > discovery_mjd + params["dt"]):
            continue

        ipix = hp.ang2pix(nside, theta=ra, phi=dec, lonlat=True)
        cumprob = prob_data_cumsum[ipix]

        transients_struct["data"] = np.append(transients_struct["data"],np.array([[ra,dec,mag,discovery_mjd,nondetection_mjd,z,mpc,cumprob]]),axis=0)
        transients_struct["name"].append(name)
        transients_struct["classification"].append(classification)

    transients_struct["name"] = np.array(transients_struct["name"])
    transients_struct["classification"] = np.array(transients_struct["classification"])

    idx = np.argsort(transients_struct["data"][:,7])
    transients_struct["data"] = transients_struct["data"][idx,:]
    transients_struct["name"] = transients_struct["name"][idx]
    transients_struct["classification"] = transients_struct["classification"][idx]

    return transients_struct


