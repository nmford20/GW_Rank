from astropy.coordinates import SkyCoord
import pandas as pd
import astropy.units as u

def Gcoordinates(distmax,distmin):

    #NOTE:need to put directory path to catalog here
    #AND if you don't already have it, need to download the GLADE+ catalog here: http://glade.elte.hu/
    #then implement a series of cuts on it, as demonstrated in this code:
    '''
    names_indx=list(range(1,40))
    names = ['col'+str(x) for x in names_indx]
    dataframe = pd.read_csv('../HOGWARTs/GLADE+.txt',dtype='str',delimiter=' ',names=names,usecols=['col4','col9','col10','col33','col11'])
    rslt = dataframe[~dataframe['col11'].isnull()]
    rslt.to_csv('../HOGWARTs/GLADE+_cut.txt',sep=' ',na_rep='NaN',header=False,index=False)
    '''
    
    #data columns are: HyperLEDA name, ra, dec, bmag, luminosity distance
    data = pd.read_csv('../HOGWARTs/GLADE+_cut.txt',sep=' ',names=['HyperLEDA','RA','Dec','Bmag','dist'], dtype={'HyperLEDA': object})
    msk1=data[['dist']]<=distmax
    msk2=data[['dist']]>=distmin
    msk3=data[['dist']]>0
    msk4=data[['dist']]!='NaN'
    msk5=data[['Bmag']]!='NaN' 
    msk6=data[['Bmag']]!='NaN' 
    msk7=data[['Bmag']]>0 
    msk=pd.concat((msk1,msk2,msk3,msk4,msk5,msk6,msk7),axis=1)
    slct=msk.all(axis=1)
    data=data.loc[slct]

    coordinates=SkyCoord(data['RA'].values*u.deg, data['Dec'].values*u.deg,data['dist'].values*u.Mpc)
    return coordinates,data