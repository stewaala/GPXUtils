# srtm data source: http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Eurasia/

# e.g. for file N34E071
#
# n[0][0]       returns the elevation at N34E071
# n[0][1200]    returns the elevation at N34E072
# n[1200][0]    returns the elevation at N35E071
# n[1200][1200] returns the elevation at N35E072

import xml.etree.ElementTree as ET
import math
import zipfile
import numpy as np
from matplotlib import pyplot as plt

a=6378137.0
b=6356752.314245
e2=0.00669437999019758
VOID_ELE = -32768
SRTM_FILE_LOC = '/Users/Alan/Data/SRTM/'
LINES_PER_CELL = 1201
MAX_GETAMAP_POINTS = 2000

def DeltaLat(degLat):

   phi = (degLat/180.0)*math.pi
   return 111132.954 - 559.822 * math.cos(2*phi) + 1.175 * math.cos(4*phi)

def DeltaLong(degLat):

   phi = (degLat/180.0)*math.pi
   return (math.pi*a*math.cos(phi))/(180*math.sqrt(1-e2*math.sin(phi)*math.sin(phi)))

def Dist(lat1, lon1, lat2, lon2):

   dLat = DeltaLat((lat1+lat2)/2)
   dLon = DeltaLong((lat1+lat2)/2)
   x=(lon2-lon1)*dLon
   y=(lat2-lat1)*dLat
   return math.sqrt(x*x + y*y)

def RouteDistance(sfile):

   tree = ET.parse(sfile)
   root = tree.getroot()
   pts = ParseRoute(root)

   SumDist = 0.0
   for i in range(1,len(pts)):
      lat1=float(pts[i-1]['lat'])
      lon1=float(pts[i-1]['lon'])
      lat2=float(pts[i]['lat'])
      lon2=float(pts[i]['lon'])
      SumDist+=Dist(lat1,lon1,lat2,lon2)/1000.0

   return SumDist
        
def NodeType(sTag):

   i = sTag.find('}')
   return sTag[i+1:].upper()

def ParseRoute(node):

   s=NodeType(node.tag)

   if s in ('GPX', 'TRK', 'TRKSEG', 'RTE', 'RTESEG'):

      pts = []
      for child in node:
         x=ParseRoute(child)
         if type(x) is dict:
            pts.append(x)
         elif type(x) is list:
            pts += x
      return pts

   elif s in ('TRKPT', 'RTEPT'):

      d={}
      d['lat']=float(node.attrib['lat'])
      d['lon']=float(node.attrib['lon'])
      for c in node:
         if NodeType(c.tag)=='ELE':
            d['ele']=float(c.text)
            break
      return d
   
   else:
      return None

def GetAscent(sfile):

   tree = ET.parse(sfile)
   root = tree.getroot()
   pts = ParseRoute(root)

   ascent=0.0
   for x2, x1 in zip(pts[1:], pts[:-1]):
      if 'ele' in x1 and 'ele' in x2:
         ascent += max(x2['ele'] - x1['ele'],0)
   return ascent

def ChartProfile(sfile):

   # the bit using srtm elevation doesn't work for some reason - the bit using the native ele in the file does not work

   tree = ET.parse(sfile)
   root = tree.getroot()
   pts = ParseRoute(root)

   n=len(pts)
   x=[0]*n
   y=[0]*n
   ncoords = np.arange(n*2).reshape(n,2)

   cdist=0
   for i in range(n):

      if i==0:
         cdist=0
      else:
         cdist+=Dist(pts[i-1]['lat'], pts[i-1]['lon'], pts[i]['lat'], pts[i]['lat'])

      x[i]=cdist
      y[i]=pts[i]['ele']

      ncoords[i,0]=float(pts[i]['lat'])
      ncoords[i,1]=float(pts[i]['lon'])

   print ncoords
   srtmelev=GetListElev(ncoords)
   
   plt.clf()
   plt.plot(x,srtmelev)
   plt.show()
   

def SRTM_File(lat, lon):

   return ('N' if lat >= 0 else 'S') + "{0:02n}".format(abs(math.floor(lat))) + ('E' if lon>=0 else 'W') + "{0:03n}".format(abs(math.floor(lon)))

def LoadTile(s):

   sfile = SRTM_FILE_LOC + s + '.hgt.zip'
   zip_file = zipfile.ZipFile(sfile,'r')
   zip_file_name = zip_file.namelist()[0]
   hgt_string = zip_file.read(zip_file_name)
   zip_file.close()

   return np.flipud(((np.fromstring(string=hgt_string, dtype='int16')).byteswap()).reshape(LINES_PER_CELL,LINES_PER_CELL))

def GetRowCol(lat, lon):

   x=lat-int(lat)
   if x<0: x=1+x
   nrow=round(x*LINES_PER_CELL-1)

   x=lon-int(lon)
   if x<0: x=1+x
   ncol=round(x*LINES_PER_CELL-1)

   return (nrow, ncol)

def GetElev(lat, lon):

   s=SRTM_File(lat, lon)
   nTile = LoadTile(s)
   (nrow, ncol) = GetRowCol(lat, lon)
   
   return nTile[nrow][ncol]

def GetListElev(npCoord):

   res = [0]*npCoord.shape[0]
   d={}
   i=0

   for lat, lon in npCoord:

      s=SRTM_File(lat, lon)

      if s in d:
         nTile = d[s]
      else:
         nTile = LoadTile(s)
         d[s]=nTile

      (nrow, ncol) = GetRowCol(lat, lon)

      x=nTile[nrow][ncol]
      res[i] = None if x==VOID_ELE else x
      i+=1

   return res

def ElevDemo(lat, lat_from, lat_to):
   
   # plot elevation profile at constant latitude=lat from lat_from to lat_to degrees longitude

   n=200
   x=np.array([[lat]*n,np.linspace(lat_from,lat_to,n)]).transpose()

   y=GetListElev(x)
   plt.clf()
   plt.plot(x[:,1],y)
   plt.show()

def CleanGPXFile(sfile):

   tree = ET.parse(sfile)
   root = tree.getroot()
   pts = ParseRoute(root)

   filt = 1+ len(pts)/MAX_GETAMAP_POINTS

   outfile = sfile[:-4] + '_' + str(filt) + '.gpx'
   f = open(outfile, 'w')
   f.writelines('<?xml version = "1.0" encoding="UTF-8" ?>\n')
   f.writelines('<gpx xmlns="http://www.topografix.com/GPX/1/1" creator = "alan.stewart" version = "1.1">\n')
   f.writelines('<trk><trkseg>\n')
   for i, x in enumerate(pts):
      if i % filt == 0:
         s = '   <trkpt lat="' + "%.8f" % x['lat'] + '" lon="' + "%.8f" % x['lon'] + '">'
         if 'ele' in x: s+= '<ele>' + "%.2f" % x['ele'] + '</ele>'
         s += '</trkpt>'
         f.writelines(s + '\n')

   f.writelines('</trkseg></trk>\n')
   f.writelines('</gpx>\n')
   f.close()


        
                        
