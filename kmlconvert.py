# Version 20150128
# D Burk
# Michigan State
#
# KMLconvert brings in the MSU HYPOUT.ASC data format file for hypocenter relocations and creates a Google Earch
# KML file of markers for all of the epicenters. It's still a crude bit of code that's designed for a special
# purpose for a very specialized data format.

# Now, the most important part -- The legalese:
# COPYRIGHT ©  BOARD OF TRUSTEES OF MICHIGAN STATE UNIVERSITY
# ALL RIGHTS RESERVED

# PERMISSION IS GRANTED TO USE, COPY, COMBINE AND/OR MERGE, CREATE DERIVATIVE
# WORKS AND REDISTRIBUTE THIS SOFTWARE AND SUCH DERIVATIVE WORKS FOR ANY PURPOSE,
# SO LONG AS THE NAME OF MICHIGAN STATE UNIVERSITY IS NOT USED IN ANY ADVERTISING
# OR PUBLICITY PERTAINING TO THE USE OR DISTRIBUTION OF THIS SOFTWARE WITHOUT 
# SPECIFIC, WRITTEN PRIOR AUTHORIZATION.  IF THE ABOVE COPYRIGHT NOTICE OR ANY
# OTHER IDENTIFICATION OF MICHIGAN STATE UNIVERSITY IS INCLUDED IN ANY COPY OF 
# ANY PORTION OF THIS SOFTWARE, THEN THE DISCLAIMER BELOW MUST ALSO BE INCLUDED.

# THIS SOFTWARE IS PROVIDED AS IS, WITHOUT REPRESENTATION FROM MICHIGAN STATE
# UNIVERSITY AS TO ITS FITNESS FOR ANY PURPOSE, AND WITHOUT WARRANTY BY MICHIGAN
# STATE UNIVERSITY OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT
# LIMITATION THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE.

# THE MICHIGAN STATE UNIVERSITY BOARD OF TRUSTEES SHALL NOT BE LIABLE FOR ANY
# DAMAGES, INCLUDING SPECIAL, INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
# WITH RESPECT TO ANY CLAIM ARISING OUT OF OR IN CONNECTION WITH THE USE OF
# THE SOFTWARE, EVEN IF IT HAS BEEN OR IS HEREAFTER ADVISED OF THE POSSIBILITY
# OF SUCH DAMAGES.

import sys
import re
import simplekml
import numpy as np
import math
 
def dist_on_unit_sphere(lat1, long1, lat2, long2):
 
    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
         
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
         
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
         
    # Compute spherical distance from spherical coordinates.
         
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta, phi)
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
     
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )
 
    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    return arc
 

# Arguments include the file name of the events file and the target file name.

kml = simplekml.Kml()
style1 = simplekml.Style()
style1.labelstyle.color = simplekml.Color.blue  # Make the text green
style1.labelstyle.scale = 0.7  # Make
style1.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
style2 = simplekml.Style()
style2.labelstyle.color = simplekml.Color.yellow  # Make the text red
style2.labelstyle.scale = 0.7  # Make
style2.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
style3 = simplekml.Style()
style3.labelstyle.color = simplekml.Color.red  # Make the text red
style3.labelstyle.scale = 0.7  # Make
style3.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
style4 = simplekml.Style()
style4.labelstyle.color = simplekml.Color.green  # Make the text red
style4.labelstyle.scale = 1.0  # Make
style4.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/target.png'


uncert_btm = 1.0

if (len(sys.argv) == 3):
    fout = sys.argv[2]
    file = sys.argv[1]
    flist = sys.argv[2]+".txt"
elif (len(sys.argv) == 2):
  fout = "c://seismo/output.kml"
  flist = "c://seismo/output.txt"
  file = sys.argv[1]    
elif (len(sys.argv) ==1):      
  file = "c://EVENTCAT//Locations//hypout2.asc"
  fout = "c://EVENTCAT//locations//hypocent2.kml"
  flist = "c://EVENTCAT//locations//hypocent2.txt"
else:
  print "useage: kmlconvert EVENTFILE TARGETFILE"
  sys.exit()

print file," outputs to ",fout

infile = open(file,'r').readlines()
input = []
depths = []
uncerts = []
resids = []
latitudes = []
longitudes = []
locdev = []     # This is used to tally up distances from the average centerpoint of all the events in file

hypocenters = open("c://seismo/hypocent","w")

#                        Parse out the fifth iteration of the hypocenter location
for line in infile:
  if (line[0] == '5'):
    hh = line[3:5]
    mm = line[6:8]
    ss = line[9:14]
    lat = line[23:30]
    lon = line[41:50]
    depth = line[62:67]
    depth_uncert = line[68:73]
    residual = line[75:80]

    
    if (residual =="*****"):
      residual = "99.99"
   
#                       Convert it to KML if it passes the test. 
#                       Depth_uncert must be less than 10.
#                       Residuals should be less than 1
    uncert = float(depth_uncert.strip(' "'))
    resid = float(residual.strip(' "'))
    flat = float(lat.strip(' "'))
    flon = float(lon.strip(' "'))



    if (uncert < 100.0) and (resid < 10.0):
      depths.append(float(depth))
      uncerts.append(uncert)
      resids.append(resid)
      latitudes.append(flat)
      longitudes.append(flon)

#      print "%s:%s:%s lat: %s Lon: %s Depth: %s Uncert: %s Resid: %s" % (hh,mm,ss,lat,lon,depth,depth_uncert,residual)
      pnt = kml.newpoint(name=depth, coords=[(flon,flat)])
      pnt.lookat = simplekml.LookAt(gxaltitudemode=simplekml.GxAltitudeMode.relativetoseafloor, latitude = flat, longitude = flon, range=1000, heading=360, tilt=0.1)
      if (float(depth)< 0.1 ):
        pnt.style = style2 # yellow
      elif (float(depth)>15.0 ):
        pnt.style = style3 # red
      else:
        pnt.style = style1
    else:
        print "****%s:%s:%s lat: %s Lon: %s Depth: %s Uncert: %s Resid: %s" % (hh,mm,ss,lat,lon,depth,depth_uncert,residual) 
# Create the analysis output.
latavg = np.mean(latitudes)
lonavg = np.mean(longitudes)
cluster = "Cluster avg depth = "+str(np.mean(depths))
pnt = kml.newpoint(name=cluster, coords=[(lonavg,latavg)])
pnt.lookat = simplekml.LookAt(gxaltitudemode=simplekml.GxAltitudeMode.relativetoseafloor, latitude = flat, longitude = flon, range=1000, heading=360, tilt=0.1)
pnt.style = style4
kml.save(fout)
      # Calculate the distances from each epicenter to the mean epicenter
      # Then tally up the distances and divide by the number of points.
      # As the epicenters converge, the value should decrease.
clusterdev = []
for i in range(0,len(latitudes)):  # find the distance from the mean epicenter to each epicenter
    clusterdev.append(dist_on_unit_sphere(latitudes[i], longitudes[i], latavg, lonavg))
avgdistance = 6373*np.sum(clusterdev)/len(clusterdev)

print "The average distance of events to mean epicenter is {0:.3f} km ".format(avgdistance)
# print clusterdev
# count the number of depths equal to zero
    
print "0 km depths: {0}.    33km depths: {1}    mean depth:  {2:2f} Km \n".format(depths.count(0.0),depths.count(33.0),np.mean(depths))
print "Mean Residual: {0:2f}        Maximum residual: {1:2f}".format(np.mean(resids),np.max(resids))
print "Mean uncertainty: {0:2f}     Max uncert: {1:2f}".format(np.mean(uncerts),np.max(uncerts))


# infile.close()
# output.close()

    






