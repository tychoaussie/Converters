__author__ = "Daniel Burk <burkdani@msu.edu>"
__version__ = "20150811"
__license__ = "MIT"

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


# Interim version - works only on PARxCH files at this time. - drb

# Version 2: Convert all DAT files, regardless of them being from the USBXCH or the PARXCH
# digitizer, regardless if they are DAT, ASC, ASC FW, or CSV format.
#
# Fixed bug in convert definition that didn't like being in the same working directory as the data.
#
# Partially written - Need to create the par1 parser, usb3 parser, and usb4 parser before operation.
#



import os
import sys
import numpy as np
import csv
import time                                # we use the sleep function to enable conversion of the DAT
import string                              # We use the string.find utility for file naming
import datetime                            # Used for assembling the time stamps in PARxCH ascii files
import subprocess                          # used for the execution of command-line programs
from obspy.core import read, Trace, Stream, UTCDateTime
from obspy.sac import SacIO


class dat2sac(object):
    '''par2sac.py is a utility for batch converting Symmetric Research DAT files into sac files
       for a whole directory. 

       The program will convert non-sequential DAT files into separate csv files, and all
       sequential files will be concatenated into a single csv file. It is therefore important
       to rename or move any sequential csv file that you do not wish to include within the csv.
        
       Two critical requirements are :
       dat2asc.exe (2010 build) exists in the directory: "c:\Python27\dat2asc.exe"
       par42asc.exe exists in the directory: "c:\Python27\par42asc.exe"

       Both of these are available in the GITHUB repository:
       http://www.github.com/tychoaussie/Converters


       Syntax: dat2csv target_directory         

       Typical useage:
       <ObsPy> C:\Python27\scripts> python Dat2sac.py c:/calibration/station/ 

    '''

def default():
    print "\n\n par2sac.py is a utility for batch converting Symmetric Research DAT files \n \
into sac files for a whole directory.\n\n \
The program will convert non-sequential DAT files into separate ASC files, and \
then convert them into standard SAC files. \n\n \
Two critical requirements are :\n \
  - dat2asc.exe (2010 build) exists in the directory: c:\Python27\dat2asc.exe\n \
  - par42asc.exe exists in the directory:             c:\Python27\par42asc.exe\n\n \
Both of these are available in the GITHUB repository:\n \
http://www.github.com/tychoaussie/Converters\n\n\n \
Syntax: dat2csv target_directory\n\n \
Typical useage:\n \
<ObsPy> C:\Python27\scripts> python Dat2sac.py c:/calibration/station/ \n"
    quit()
    

                            # Function getcal: Retrieve the calibration control file
def getcal(calcontrol): 
    # calcontrol needs to include calibration constants as well as the station name and the channel name, and the particular constant for that channel.
    # Thus a third line is necessary that specifies the channel identifier and the channel assignment of that channel.
    # Channel name is located in the top row already, and it's position is associated with the sensitivity. So the third row designates the UUT and the
    # laser position channel.
    cconstant = ["MSU","CH0",1.0,"CH1",1.0,"CH2",1.0,"CH3",1.0,1.0,1.0,0.707,1.0,1,2]
    try:
        with open(calcontrol,'r') as fin:
            list = csv.reader(fin)
            rowcnt=0            
            stack = []
            header = []
            calconstants = []
            selection = []
            for row in list:
                stack.append(row)
            header = stack[0]
            calconstants = stack[1]
            if len(stack)==3:# old calibration file format, so prompt the user for the appropriate channel assignments.
                selection = stack[2]
            else:
                for i in range(0,4):
                    print "\n\nNo channels selected for calibration in cal control file."
                for i in range(0,4):
                    print"Note: Channel {0} is listed as channel number {1}".format(header[i+1],i)
            
                selection.append(int(raw_input('Choose the channel number representing the unit under test:  ')))
                selection.append(int(raw_input('Choose the channel number representing the laser position sensor:  ')))

            cconstant[0] = calconstants[0]        # (test) Station name
            cconstant[1] = header[1]              # (text) Channel name for CH0
            cconstant[2] = float(calconstants[1]) # (float) adccal[0]: cal constant for ch 0 (microvolts / count)
            cconstant[3] = header[2]              # (text) Channel name for CH1
            cconstant[4] = float(calconstants[2]) # (float) adccal[1]: cal constant for ch 1 (microvolts / count)
            cconstant[5] = header[3]              # (text) Channel name for CH2
            cconstant[6] = float(calconstants[3]) # (float) adccal[2]: cal constant for ch 2 (microvolts / count)
            cconstant[7] = header[4]              # (text) Channel name for CH3
            cconstant[8] = float(calconstants[4]) # (float) adccal[3]: cal constant for ch 3 (microvolts / count)
            cconstant[9] = float(calconstants[5]) # (float) laserres: cal constant for the laser ( mV / micron)
            cconstant[10] = float(calconstants[6])# (float) lcalconst: cal constant for geometry correction factor
            cconstant[11] = float(calconstants[7])# (float) h: Damping ratio for the seismometer as measured by engineer.
            cconstant[12] = float(calconstants[8])# (float) resfreq: Free period resonance freq. as measured by engineer.
            cconstant[13] = int(selection[0])     # channel number of channel being tested
            cconstant[14] = int(selection[1])     # channel number of the laser position sensor data
    except:
        print "Calibration control file not found!\n"
        print "You may enter the calibration information manually, or else press <ctrl>C"
        print "then place the cal control file in the directory containing the data files."
        cconstant = ["","",1.0,"",1.0,"",1.0,"",1.0,1.0,1.0,1.0,1.0,2,3]
        try:
            cconstant[0] = raw_input('Enter the station name: ')# (test) Station name
            print "Enter the four channel names, starting with channel 0:\n"
            cconstant[1] = raw_input('Channel 0: ') # (text) Channel name for CH0
            cconstant[3] = raw_input('Channel 1: ') # (text) Channel name for CH1
            cconstant[5] = raw_input('Channel 2: ') # (float) adccal[1]: cal constant for ch 1 (microvolts / count)
            cconstant[7] = raw_input('Channel 3: ') # (text) Channel name for CH2
            print "Enter the channel calibration constants for the above four channels.\n"
            cconstant[2] = raw_input('CH0 calibration constant (uV/count): ') # (float) adccal[2]: cal constant for ch 2 (microvolts / count)
            cconstant[4] = raw_input('CH1 calibration constant (uV/count): ')              # (text) Channel name for CH3
            cconstant[6] = raw_input('CH2 calibration constant (uV/count): ') # (float) adccal[3]: cal constant for ch 3 (microvolts / count)
            cconstant[8] = raw_input('CH3 calibration constant (uV/count): ')  # (float) laserres: cal constant for the laser ( mV / micron)
            cconstant[9] = raw_input('Enter the resolution for the laser position sensor (mV/micron): ')
#            print "The laser geometry correction ratio is the ratio of distance from pendulum pivot to center of mass"
#            print "divided by the distance from pendulum pivot to the measurement point of the laser beam."
            cconstant[10] = 1.0 #raw_input('What is the laser geometry correction constant ratio? ') # (float) lcalconst: cal constant for geometry correction factor
            cconstant[11] = 0.707 #raw_input('Enter the measured damping ratio for the seismometer under test. ') # (float) h: Damping ratio for the seismometer as measured by engineer.
            cconstant[12] = 1.0 #raw_input('Enter the resonance frequency of the seismometer under test. (in Hz). ') # (float) resfreq: Free period resonance freq. as measured by engineer.
            cconstant[13] = 1 #raw_input('Enter the channel number (0 through 3) for the seismometer under test. ')    # channel number of channel being tested
            cconstant[14] = 2 #raw_input('Enter the channel number (0 through 3) repesenting the laser position sensor. ')     # channel number of the laser position sensor
        except:
           print "Error during manual input of the calibration constant parameters.\n"
           print "Setting the paramters to default settings."
           cconstant = ["MSU","CH0",1.0,"CH1",1.0,"CH2",1.0,"CH3",1.0,1.0,1.0,0.707,1.0,1,2]
    return(cconstant)






def location(gpgga):
#                                   Compute the mean gps location from the data
    lat = []
    lon = []
    for i in range(0,len(gpgga)):
        latp = 0
        lonp = 0
        if gpgga[i][3] == 'N':
            latp = 1
        else:
            latp = -1
        if gpgga[i][5] == 'E':
            lonp = 1
        else:
            lonp = -1
        flat = float(gpgga[i][2][0:2])+float(gpgga[i][2][2:10])/60
        flon = float(gpgga[i][4][0:3])+float(gpgga[i][4][3:11])/60
        lat.append(flat*latp)
        lon.append(flon*lonp)
    latitude = np.mean(lat)
    longitude = np.mean(lon)
    return(latitude,longitude)





def par1(stack):
#                             Build the par1 stack parsing routine here for USB4CH .csv formatted data

    return(header,data,stinfo)





def par2(stack): # File format consists of about 116 lines of header, then data that is interspersed with GPS information.
                  # Header is 52 lines of misc information, followed by four "time point" records.
                  # Line 0 is the file origin
                  # Header begins with the word "Header" at line 2
                  # Information important for conversion found within the header:
                    # 'ChannelsAnalog' = 4 (typically)
                    # 'ChannelPtsPerRecord' = 128 lists the number of samples we'll find before breaking for some GPS info
                    # 'ChannelPtsPerFile' = 5120 lists the number of samples to be found within the file
                    # 'Channel 00' = CH00 Channel name for channel 0
                    # 'Channel 01'                         channel 1
                    # 'Channel 02'
                    # 'Channel 03'
                    # 'Channel 04' = MARK = The GPS channel that marks the 0.00 second sample for each and every second.
                    # 'AtodSamplingRate' = 30.001921 lists the assumed sample rate of the digitizer
                  # TimePoint[0]:   There are four time points listed in the header. The essential information is:
                    # 'SampleFromRun' = The corresponding sample ID that coincides with the following time point.
                    # 'Year'
                    # 'Month'
                    # 'Day'
                    # 'Hour'
                    # 'Minute'
                    # 'Second'
                    # 'Microsecond'
                  #   Timepoint will be used to time-code the first sample within the file.
                    # After time-coding the sample, we find the sample ID within the data series corresponding to this sample.
                    # This should be the first sample listed within the file. However if it is not, we can back-calculate time
                    # for preceding samples using the established sample rate.
                    # We can also establish file integrity by looking at the 'mark' to make sure it's non-sero on the even
                    # second samples.
                  #   NMEA GPS data contains the station coordinates on the line tagged GPGGA: It is comma separated.
                    # $GPGGA,194614,4435.7663,N,03327.4974,E,1,03,6.9,,M,30.6,M,,*63
                    # Means time taken at 19:46:14 gmt, coordinates 44 deg, 35.7663 " North, 33 deg, 27.4974" East.
                    # This will be saved to stla and stlo within the sac header.
                    
    header = []                
    data = []
    databuff = []
    exceptions = 0
    gpgga = []
    stinfo = []
    for i in range(0,len(stack)):  # Split out the data from the header and gps info
        try:
            if len(stack[i]) >1:
                # check for GPGGA strings
                if '$GPGGA' in stack[i][0]:
                    gpgga.append(stack[i])
            else:    
                line = str.split(stack[i][0])  # Break the single-element list at element i into a multi-element list
                if len(line)>=7:
                    if "Total" not in stack[i][0]: # toss out any data lists that consist only of header information
                        databuff = []
                        for j in range(2,7):
                            databuff.append(int(line[j]))
                        data.append(databuff)
                else:
                    header.append(line) # place everything else into the header, which has been cleaned of the data.
        except:
            # invalid element in the stack; likely a blank space.
            #print "warning: blankspace at line {}".format(i)
            exceptions = exceptions +1 # Just count the exceptions for now and keep trying until the end of the list.

            
    latitude,longitude = location(gpgga)  # Get an averaged location from the gpgga data.


#    print "Mean location from GPS tags: {} , {}\n".format(latitude,longitude)
    stinfo.append('MI')     # 0 default network name
    stinfo.append('MSU')    # 1 Station name
    stinfo.append(latitude) # 2 Station latitude
    stinfo.append(longitude)# 3 Station longitude
    stinfo.append(header[25][3]) # 4 Channel 0 name
    stinfo.append(header[26][3]) # 5 
    stinfo.append(header[27][3]) # 6
    stinfo.append(header[28][3]) # 7
    stinfo.append(header[29][3]) # 8 Channel 4 name
    stinfo.append(1/np.float(header[33][2])) # 9 delta is the iverse of sample rate
    stinfo.append(header[24][2]) # 10 Number of points per file
    stinfo.append(header[56][2]) # 11 Year
    stinfo.append(header[57][2]) # 12 Month
    stinfo.append(header[58][2]) # 13 Day
    stinfo.append(header[60][2]) # 14 Hour
    stinfo.append(header[61][2]) # 15 Min
    stinfo.append(header[62][2]) # 16 Sec
    stinfo.append(header[63][2]) # 17 Microsec
    
    return(header,data,stinfo)




def usb3(stack):
#                             Build the usb3 stack parsing routine here for USB4CH .csv formatted data
#                             This definition is NOT finished: Simulate it in ipython notebooks to extract times
#                             and setup the stinfo, data, and header information from a csv file.

    header = []                
    data = []
    exceptions = 0
    gpgga = []
    stinfo = []
    Channel = ["","","",""]

    units = ['Counts  ','Counts  ','Counts  ','Counts  ']
    comment = ['Velocity','Velocity','Velocity','Velocity']
    header = stack[0]

    datetime = stack[1][13]+","+stack[1][14]
    Frac_second = float(stack[1][15])
    St_time = time.strptime(datetime,"%Y/%m/%d,%H:%M:%S")

    cconstant = getcal(calcontrol) # This will need to be fixed when using targeted directories.

    Station = cconstant[0]
    outfile = infile[0:string.rfind(infile,'.')]
    for i in range(0,4):
        Channel[i]=cconstant[(2*i)+1]

    Samplecount = len(stack)
    print "Sample count stands at {} samples.".format(Samplecount)
    Delta = (float(stack[len(stack)-1][12])-float(stack[1][12]))/len(stack)
    sacfile = outfile[:string.find(infile,'.')]+'{}'.format(i)+'.sac'
        #
        # stack[1] = channel 1 time history
        # .
        #
        # stack[4] = channel 4 time history
        #
    stinfo.append('MI')     # 0 default network name
    stinfo.append('MSU')    # 1 Station name
    stinfo.append(latitude) # 2 Station latitude
    stinfo.append(longitude)# 3 Station longitude
    stinfo.append(header[25][3]) # 4 Channel 0 name
    stinfo.append(header[26][3]) # 5 
    stinfo.append(header[27][3]) # 6
    stinfo.append(header[28][3]) # 7
    stinfo.append(header[29][3]) # 8 Channel 4 name
    stinfo.append(1/np.float(header[33][2])) # 9 delta is the iverse of sample rate
    stinfo.append(header[24][2]) # 10 Number of points per file
    stinfo.append(header[56][2]) # 11 Year
    stinfo.append(header[57][2]) # 12 Month
    stinfo.append(header[58][2]) # 13 Day
    stinfo.append(header[60][2]) # 14 Hour
    stinfo.append(header[61][2]) # 15 Min
    stinfo.append(header[62][2]) # 16 Sec
    stinfo.append(header[63][2]) # 17 Microsec
    return(header,data,stinfo)



def usb4(stack):
#                             Build the usb3 stack parsing routine here for USB4CH .csv formatted data

    return(header,data,stinfo)


 

# def load(infile): # Load the csv file
#    with open(infile,'r') as fin:
#        list = csv.reader(fin)
#        rowcnt=0
#        stack = []
#        header = []
#        for row in list:
#
#          Bring in the data and create a list of lists, each of which
#           corresponds with a given sample.
#
#            if rowcnt == 0:
#                header.append(row)
#            else:
#                stack.append(row)
#            rowcnt = 1
#    return(header,stack)



def load(infile): # Load the csv file
    # Identify the file type by looking at specific fields. This is a Symmetric Research converter.
    # Filetype 1: PAR4CH with no timing associated with the sample
    # Filetype 2: PAR4CH with timing associated with the samples
    # Filetype 3: USB4CH in CSV format
    
    with open(infile,'r') as fin: # load the data into a stack to determine file type.
        list = csv.reader(fin)
        rowcnt=0
        ftype = 0   # Unknown data format.
        stack = []
        header = []
        data = []
        gpgga = []
        for row in list:
            stack.append(row)
        if ('PAR4CH' in stack[7][0]):
            ftype = 1     # This file doesn't include timing for each sample and is a PAR4CH data format
            if(len(stack[119][0]) > 80):
                ftype = 2 # This file includes sample timing and is also a PAR4CH data format
                header,data,stinfo = par2(stack)
            else:
                header,data,stinfo = par1(stack)
        elif ('Pt#,' in stack[0][0]): # The file is a Symres USB4CH csv data file format
            ftype = 3
            header,data,stinfo = usb3(stack)
        elif ('Dat2Asc' in stack[0][0]):
            ftype = 4      # This is a fixed width Symres USB4CH data file format ascii file
            header,data,stinfo = usb4(stack)
        else:
            print "Unknown data format."
            data = stack
        
            
    print "The file type = {}".format(ftype)
    return(ftype,header,data,stinfo,stack)



                        # process dat2sac
                        # Converts a dat file into a sac file format

def asc2sac(infile,netnames):
#                              These are defaults to be overwritten in either cconstant or stinfo.
# stinfo[0] = Network name ( prompt from user)
# stinfo[1] = station name ( prompt from user)
# stinfo[2] = station latitude ( taken from nmea data)
# stinfo[3] = station longitude ( taken from nmea data)
# stinfo[4] = channel 0 name ( found in header at header[26][3]) 
# stinfo[5] = channel 1 name header[27]
# stinfo[6] = channel 2 name header[28]
# stinfo[7] = channel 3 name
# stinfo[8] = channel 5 name (mark)
# stinfo[9] = delta, which is the sample period
# stinfo[10] = number of samples
# stinfo[11] = start year
# stinfo[12] = start day
# stinfo[13] = start hour
# stinfo[14] = start minute
# stinfo[15] = start second
# stinfo[16] = start fractions of second
 
    Channel = ["","","","",""]
    units = ['Counts  ','Counts  ','Counts  ','Counts  ','Counts  ']
    comment = ['Velocity','Velocity','Velocity','Velocity','microvlt']
    idep = [4,4,4,4,1]


#    cconstant = getcal(calcontrol)


    (ftype,header,data,stinfo,stack) = load(infile) # Load function will search out the type of ascii file.  
#    time is now handled in stinfo and is managed in the load definition.
#    data contains the list of time history for the channels 
    nzyear = int(stinfo[11])
    nzjday = int(datetime.datetime.strptime(stinfo[11]+stinfo[12]+stinfo[13], '%Y%m%d').timetuple().tm_yday)
    nzhour = int(stinfo[14])
    nzmin = int(stinfo[15])
    nzsec = int(stinfo[16])
    nzmsec = int(int(stinfo[17])/1000)
    Delta = stinfo[9]
    Samplecount = stinfo[10]
    Network = netnames[0]
    Station = netnames[1]

    for i in range(0,len(data[0])):        # Channel names are now taken from stinfo in the load definition
        Channel[i]=stinfo[4+i]

    outfile = infile[0:string.rfind(infile,'.')]



    print "Sample count stands at {} samples.".format(Samplecount)

    sacfile = outfile[:string.find(infile,'.')]+'{}'.format(i)+'.sac'
        #
        # stack[1] = channel 1 time history
        # .
        #
        # stack[4] = channel 4 time history
        #
        #    cconstant[0] = calconstants[0]        # (test) Station name
        #    cconstant[1] = header[1]              # (text) Channel name for CH0
        #    cconstant[2] = float(calconstants[1]) # (float) adccal[0]: cal constant for ch 0 (microvolts / count)
        #    cconstant[3] = header[2]              # (text) Channel name for CH1
        #    cconstant[4] = float(calconstants[2]) # (float) adccal[1]: cal constant for ch 1 (microvolts / count)
        #    cconstant[5] = header[3]              # (text) Channel name for CH2
        #    cconstant[6] = float(calconstants[3]) # (float) adccal[2]: cal constant for ch 2 (microvolts / count)
        #    cconstant[7] = header[4]              # (text) Channel name for CH3
        #    cconstant[8] = float(calconstants[4]) # (float) adccal[3]: cal constant for ch 3 (microvolts / count)
        #    cconstant[9] = float(calconstants[5]) # (float) laserres: cal constant for the laser ( mV / micron)
        #    cconstant[10] = float(calconstants[6])# (float) lcalconst: cal constant for geometry correction factor
        #    cconstant[11] = float(calconstants[7])# (float) h: Damping ratio for the seismometer as measured by engineer.
        #    cconstant[12] = float(calconstants[8])# (float) resfreq: Free period resonance freq. as measured by engineer.
        #    cconstant[13] = int(selection[0])     # channel number of channel being tested
        #    cconstant[14] = int(selection[1])     # channel number of the laser position sensor data

    for i in range(0,len(data[0])): # Build each channel
        t = SacIO()
        b = np.arange(len(data),dtype=np.float32)   #   Establishes the size of the datastream
        for n in range(len(data)):        #   Load the array with time-history data
            b[n] = np.float32(data[n][i]) #   Convert the measurement from counts to volts.
        t.fromarray(b)

             #                     set the SAC header values
        t.SetHvalue('scale',1.00) # Set the scale for each channel. This one is important to declare.
        t.SetHvalue('delta', Delta)
        t.SetHvalue('nzyear',nzyear)
        t.SetHvalue('nzjday',nzjday)
        t.SetHvalue('nzhour',nzhour)
        t.SetHvalue('nzmin',nzmin)
        t.SetHvalue('nzsec', nzsec)
        t.SetHvalue('nzmsec', nzmsec)
        t.SetHvalue('kstnm',Station)
        t.SetHvalue('kcmpnm',Channel[i])
        t.SetHvalue('idep',idep[i]) # 4 = units of velocity (in Volts)
                              # Dependent variable choices: (1)unknown, (2)displacement(nm), 
                              # (3)velocity(nm/sec), (4)velocity(volts), 
                              # (5)nm/sec/sec
        t.SetHvalue('kinst',comment[i])       # Instrument type
        t.SetHvalue('knetwk',Network)         # Network designator
        t.SetHvalue('kuser0',units[i])        # Place the system of units into the user text field 0
        t.WriteSacBinary(outfile+"_{}.sac".format(Channel[i]))
        print " File successfully written: {0}_{1}.sac".format(outfile,Channel[i])       

  


def convert(infile,digitype,netnames):
    print infile
    ext = '.asc'
    if digitype == 'csv':
        ext = '.csv'
# digitype must determine the correct file name extension. 'csv for usb4ch, .asc for par4ch
# Need to fix this subroutine.

    target = infile[string.rfind(infile,"\\")+1:string.find(infile,'.')]+ext
    outfile = infile[:string.find(infile,'.')]+ext    
    calcontrol = infile[:string.rfind(infile,"\\")+1]+"calcontrol.cal"
    if digitype == 'par':
        subprocess.call(["c:\\Python27\\par42asc.exe",infile,"/tymd"])
    else:
        dat2csvfile = infile[:string.rfind(infile,"\\")+1]+"Dat2asc-301-Data.csv"
        subprocess.call(["c:\\Python27\\dat2asc.exe",infile,"csv"])
        print "convert {} to: \n".format(dat2csvfile)
        subprocess.call(["ren",dat2csvfile,target],shell=True)
    print 'file output to: {}'.format(outfile)
    asc2sac(outfile,netnames)


def main():
                                      #           MAIN PROGRAM BODY
                                      #  Parse the command line switches
                                      # Commmand example: c:\Python27>Python.exe Sigcal.py c:\seismo\caldata\momo
                                      # where momo is the working directory containing the csv files
                                      # as well as the calibration control file, c:\seismo\caldta\calcontrol.csv

                                      
    optioncount = len(sys.argv)
    outputfile_defined = False
    filelist = []
    dir=""
    infile = ""
    directory = os.getcwd()
    directory = directory.replace("/","\\")
#    print "This is the directory name as taken from the computer:",directory
    digitype = 'par' # The default data type
    if optioncount > 1:
        for i in range(1,optioncount):
#            print 'test argument {0} = {1}'.format(i,sys.argv[i])
            if '.' not in sys.argv[i] and ("/" in sys.argv[i] or "\\" in sys.argv[i]):  # Assume CLA is a directory
                directory = sys.argv[i]
                directory = directory.replace("/","\\")               
#                print 'target directory: {}'.format(directory)
                if directory[-1:] !="\\":
                    directory = directory+"\\"
                    try:
                        filelist = os.listdir(directory)            
                    except:
                        print "Command line parameter must consist of a valid directory. {}".format(directory)
                        sys.exit(0)
                else:
                    if directory[-1:] !="\\":
                        directory = directory+"\\"
                        filelist = os.listdir(directory)
            else:
                if string.lower(sys.argv[i]) == 'par': # Is this a command to convert par4ch files?
                    digitype = 'par'
#                    print "Setting the digitype to par"
                elif string.lower(sys.argv[i]) == 'usb': # Is this a command to convert USB4CH files?
                    digitype = 'usb'
                elif '.dat' in string.lower(sys.argv[i]): # Looks like a discrete ref to a actual file
#                    print "Adding {} to the file list.".format(sys.argv[i])
                    directory = ""
                    filelist.append(sys.argv[i])
    else:
        default()

    netnames = []                                             
    netnames.append(raw_input('Enter the Network designator: '))
    netnames.append(raw_input('Enter the station code: '))
    
    if ".dat" in string.lower(filelist[0]):
        infile = directory+filelist[0]
        convert(infile,digitype,netnames)

    for n in range(1,len(filelist)):           
#        print "Parsing through the file list with element {0} = {1} ".format(n,filelist[n])                     
        if ".dat" in string.lower(filelist[n]):
            try:                    
                filenum1= int(filelist[n][string.rfind(filelist[n],'.')-8:string.rfind(filelist[n],'.')])
                try: 
                    filenum0= int(filelist[n-1][string.rfind(filelist[n-1],'.')-8:string.rfind(filelist[n-1],'.')]) # If 
                    if ((filenum1-filenum0)!=1):                    # Skip sequential files that have likely been converted with prev. file
                        infile = directory+filelist[n]
                        print "Converting: ",infile
                        convert(infile,digitype,netnames)
                except:                                             # previous file failed but this one does not.
                    infile = directory+filelist[n]     
                    print "Converting: ",infile
                    convert(infile,digitype,netnames)
            except:
                try: # one last ditch effort
                    infile = directory+filelist[n]
                    print "Converting: ",infile
                    convert(infile,digitype,netnames)
                except:
                    print "File {} does not comply to standard symmetric research naming formats and must be manually converted.".format(filelist[n])

#
# Check and run the main function here:
#
if __name__ == '__main__':
  main()
 