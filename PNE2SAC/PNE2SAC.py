# -*- coding: utf-8 -*-

'''The MIT License (MIT)



Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.'''


__author__ = "Daniel Burk <burkdani@msu.edu>"
__version__ = "20200205"
__license__ = "MIT"

# -*- coding: utf-8 -*-

# Today, on November 5th, 2014 we lost a family member to illness. We will
# always remember him.

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

import os, csv, sys, numpy as np, time, calendar, string
from obspy.io.sac import SACTrace
from obspy.core import read, Trace, Stream, UTCDateTime
from obspy.signal.detrend import polynomial
import matplotlib.pyplot as plt
from scipy import interpolate

class ASCconvertError(Exception):

    def __init__(self, error_name, error_text):
        self.args = (error_name, error_text)
        self.error_name = error_name
        self.error_text = error_text

    def __str__(self):
        return repr('{0}: {1}'.format(self.error_name, self.error_text))

class ASCconvert(object):
    '''PNE2SAC is a utility for converting text files representing MSU digitized 
       records of PNE events from the Soviet Union. 
       
       Useage: Activate Pythons ObsPy package, then from the command prompt
       specify the target file. 

       Syntax1: python PNE2SAC.py infile_name
       Syntax2: python PNE2SAC.py directoryname (file_extension)
       Syntax 2 enables you to specify a whole directory to convert multiple
       files in one go.


       Typical useage:
       <ObsPy> C:\Python27\scripts> python PNE2SAC.py C:/../infil.txt 
       <ObsPy> C:\Python27\scripts> python PNE2SAC.py C:/seismic/PNE 

    '''
# version 20200205: Move to processing the header information differently. Use the "REFTIME"
#                   to represent the relative time (on the seismogram) of the first sample in the file
#                   Then, user "STARTTIME" to indicate the relative time on the seismogram where to start
#                   the digitization process. There is no longer any need to hand-calculate the actual
#                   time within the file. PNE2SAC will apply TC at the last step to correct the relative
#                   start time into absolute start time in GMT.
#                       In addition, code is now network aware. Add the field "NETWORK" to the header of
#                       all future wavetrack output files. See below.
#                   
# version 20201028: Add the polynomial signal detrend function to remove offsets from incoming data.
#                    These offsets are likely part of the digitization process.
# version 20200122: Create a database-compliant naming convention for the SAC and miniseed output.
# version 20200111: Add clickpoint recovery and Chip waveform interpolation to the code
# along with a plot output with overlay of original waveform, clickpoints, and modified waveform.

# version 20191015: Upgrade to python 3.7, and modify the response conversion
#                   to hold amplitude to centimeters. 
# CF is set to unity in the text file. 
# Change the SAC output to displacement rather than output of velocity.

# Do the channel response removal within a dataless SEED or a SAC pole/zero conversion.

#
# Load the PNE file, extract the first six parameters which make up the header file.
# Count the remaining samples which consist of an offset time and a sample.
# 
# Calculate the true start time, and calculate the delta based on offset of sample n
# minus offset of sample zero.
# 
# Calculate the number of samples in the file.
# Return the list of header information packets and the list of data points.
#
# SAChead list consists of Delta, nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec, 
# kstnm(station name), kcmnpm(channel name).
#
# DATA consists of the list of points in centimeters that we pre-convert by applying the conversion
# factor and a multiplier of 10,000 to convert to microns of ground motion.
# 
# Two lists are returned: SAChead, and Data
#
#
# Each text file must contain a header with the following eight rows:
#
#    COMMENT Neva-2-2-Stolb-409_600dpi
#    NETWORK RY
#    STATION STB
#    COMPONENT HHZ
#    REFTIME   24_JUL_1987_02:03:00.000
#    STARTTIME 24_JUL_1987_02:05:00.000
#    CF 25000
#    TC -1.5
#
#    followed by the data fields, representing time (relative to REFTIME) and amplitude (in centimeters):
#
#    0.000000	4.360000
#    0.010000	4.359900  
#    ...etc.
#
#    Note that REFTIME is the time read from the seismogram corresponding to the first sample.
#
#    STARTTIME is the relative requested start to which we will trim the data to the output start time, and
#    typically represents the minute marker immediately preceding the onset of the first arrival.
#
#    TC is the time correction factor: A number of seconds by which we need to adjust displayed time 
#        (plus the seconds offset from first sample) to yield the true GMT time.
#    PNE2SAC now automatically adds TC to time-correct the output to true GMT. There is no longer any need
#    to manually calculate it into the start time. 
#
#    When building the header, ensure that station code is in compliance with pre-existing ISC international station code
#    and that the component match standard naming convention for IRIS/SEED channel names (i.e. SHZ,BHN,HHE, etc.)






#                                   The Defs for PNE2SAC:


def load(infile):                 # Collate the file into header information and a stack of data.
    with open(infile,'r') as fin:
        list = csv.reader(fin)
        rowcnt = 0
        stack = []
        header = []
        for row in list:
            if len(row) > 0:          # skip blank lines
                r = row[0].split()    # First seven rows are header and assigned to PNE[0]
                if rowcnt < 8:
                    header.append(r)
                    rowcnt +=1
                else:                 # Stack represents data and is assigned to PNE[1]
                    stack.append(r)
                    rowcnt+=1
                
        return(header,stack)

#
# Error during processing is usually caused by two wavetrack clickpoints on the same horizontal time slice
# It is also caused by setting the sample interval of the wavetrack too low. Best results are obtained
# when sample interval is set to 0.01 seconds in Wavetrack for a 100 sample per second output.
#

def log_error(point,x,y,error):
    if error == -1:
        print(f"Error at time = {point[0][0]} seconds where calculated slope segments are equal.")
    else:
        print(f"Error in calculation of clickpoint time: {point[0][0]:2.3f} < {error:2.3f} < {point[1][0]:2.3f} sec.")
        print(f"Time discrepency of {(error-point[0][0]):2.3f} seconds.") 
        print(f"Use Linear interpolation set sample time to {x:2.3f} seconds.\n")



#  https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines

def line_intercept(point):
                            # Points are configured as [time(x), amplitude(y), slope(m)]
    error = 0
    point1 = point[0]
    point2 = point[1]
    b1 = (point1[1] - point1[2]*point1[0])  #b = amplitude(y) - slope(m) * time(x)
    b2 = (point2[1] - point2[2]*point2[0])  # this is b: the y intercept.
    m1 = point1[2] #   
    m2 = point2[2] #
    if m1 == m2:
        error = -1
        interpolated_time = point1[0]+(point2[0]-point1[0])/2.0
        interpolated_amplitude = point1[1]+(point2[1]-point1[1])/2.0
        return (interpolated_time,interpolated_amplitude,error)
    x = (b2 - b1) / (m1 - m2) # Calculate the x position of the interection of the lines
    y = m1 * x + b1

    # If the calculated point falls more than .002 seconds outside of the wavetrack points, assume
    # there is a discontinuity within wavetrack data and interpolate.
    
    if not point1[0] <= x <= point2[0] and ((x - point1[0] < -.002) or (x - point2[0] > 0.002)) :
        error = x
        x = point1[0]+(point2[0]-point1[0])/2.0 # interpolate a midpoint
        y = point1[1]+(point2[1]-point1[1])/2.0
    return (x,y,error) # Pass the miscalculated value for error logging.




# (PCHIP stands for Piecewise Cubic Hermite Interpolating Polynomial)
# clues taken from here: https://stackoverflow.com/questions/31221444/intuitive-interpolation-between-unevenly-spaced-points

def Pchip(time,amp): # input the time and amplitude series. method describes the desired interpolation method
    # method = 0 [default] = Pchip interpolation methodology
    # method = 1 = 
    # Slope = difference in amplitude across two adjacent samples divided by time difference between the time of those two samples
    slopes = []
    bslopes = []
    fslopes = []
    deltat = [] # [time[1]-time[0]]
    for i in range(1,len(time)):
        deltat.append(time[i]-time[i-1])
    dt = np.mean(deltat) # This is the standard sample rate (IN SECONDS) from the wavetrack output.  

    # We now have time, amp, and slopes
    # Assume that the first element slope is the same as the second element's slope

    slopes.append((amp[1]-amp[0])/dt) # front load with the first point
    fslopes.append((amp[1]-amp[0])/dt) # fore-looking amplitude change/time (cm/sec)
    bslopes.append((amp[1]-amp[0])/dt) # back-looking amplitude change/time (cm/sec)

    for i in range(1,len(time)-1):
        slopes.append((amp[i]-amp[i-1])/dt) # amplitude change/time (cm/sec)
        bslopes.append((amp[i]-amp[i-1])/dt) # back-looking amplitude change/time (cm/sec)
        fslopes.append((amp[i+1]-amp[i])/dt) # fore-looking amplitude change/time (cm/sec)

    slopes.append((amp[len(amp)-1]-amp[len(amp)-2])/dt) #  load with the last point
    bslopes.append((amp[len(amp)-1]-amp[len(amp)-2])/dt) # backward-looking amplitude change/time (cm/sec)
    fslopes.append((amp[len(amp)-1]-amp[len(amp)-2])/dt) # forward-looking amplitude change/time (cm/sec)

    # We now have a full list for each point including forward and back-looking slopes for evaluation.
    # front-load list with the initial point, assumed to represent the initial mouse click.
    #
    #                    PARSE OUT THE CLICK POINTS
    #
    ds = [slopes[0]]   # delta-s, the change in slope
    dstime = [time[0]] # frontload, assuming the first point in the list is an initial pick point
    dspointer = [0]    # From the original list, this is the index for the point.
    dsamp = [amp[0]]   # frontload with the first amplitude and calculated slope: essentially, the first clickpoint.

    foreslopes = [fslopes[0]]
    backslopes = [bslopes[0]]
    points = [[[dstime[0],dsamp[0],foreslopes[0]],[]]]
    second_point = False
    for i in range(1,len(slopes)-1): # dont check the first and last item as these are automatically added to the list
        if not second_point:
            if (np.abs(slopes[i+1] - slopes[i])) > dt/2 :
                ds.append((slopes[i] - slopes[i-1]))       # from a different click point
                backslopes.append(bslopes[i])
                foreslopes.append(fslopes[i])
                dsamp.append(amp[i])
                dstime.append(time[i])
                dspointer.append(i)
                second_point = True
                point0 = [time[i],amp[i],bslopes[i]] # pointer, time, amplitude, backslope for first point
        else:
            second_point = False
            ds.append((slopes[i] - slopes[i-1]))       # ADD the second point that brackets the clickpoint
            backslopes.append(bslopes[i])
            foreslopes.append(fslopes[i])
            dsamp.append(amp[i])
            dstime.append(time[i])
            dspointer.append(i)
            point1 = [time[i],amp[i],fslopes[i]]
            points.append([point0,point1])

    ds.append(np.abs(slopes[len(slopes)-1] - slopes[len(slopes)-2]))
    dstime.append(time[len(time)-1])
    dspointer.append(len(time)-1)
    dsamp.append(amp[len(amp)-1])
    foreslopes.append(fslopes[len(fslopes)-1])
    backslopes.append(bslopes[len(bslopes)-1])        
    endpoint = len(dstime)-1
    points.append([[dstime[endpoint],dsamp[endpoint],backslopes[endpoint]],[]]) # Add the last point, (last clickpoint)

    # At this point, we have a list of points that bracket the actual click point. 
    # Now, process this list of points to yield a list of click points and the times represented by the click point.
    # Now iterpolate each point and build a list of times and amplitudes

    clicktime = []
    clickamp = []
    errtime = []
    for point in points:
        if not point[1]:       # a second point does not exist for the first and last point within the list.
            clicktime.append(point[0][0])
            clickamp.append(point[0][1])
        else:                  # find the intercept point from the points list and add it to the final list.
            x,y,error = line_intercept(point)

            clicktime.append(x)
            clickamp.append(y)
            if error != 0:
                log_error(point,x,y,error) # Process the error
                errtime.append(x) # Record whenever a linerly interpolated point exists.

    # At this point, clicktime represents the time of each click point and clickamp represents amplitude of each mouse click

    x_data = np.array(clicktime)
    y_data = np.array(clickamp)

    Pchip_time = np.linspace(min(x_data), max(x_data), len(slopes))
    bi = interpolate.PchipInterpolator(x_data, y_data)
    Pchip_amplitude = bi(Pchip_time)
    return(Pchip_time,Pchip_amplitude,clicktime,clickamp,errtime)




def main():
#                               Parse the command line switches
    optioncount = len(sys.argv)
    SAC = False
    CFflag = False
    outputfile_defined = False
    filelist = []
    dir=""
    extension = '.txt'
    if optioncount > 1:

        if optioncount == 4:
            if sys.argv[3] == '-s':
                SAC = True
            if sys.argv[3] == '-cf':
                CFflag = True
            outfile = sys.argv[2]
            infile = sys.argv[1]
            filelist.append(infile)

        elif optioncount == 3:
            if "." in sys.argv[1]:
                
                infile = sys.argv[1]
                filelist.append(infile)
                outputfile_defined = True
                outfile = sys.argv[2]
            else:
                if len(sys.argv[2])==4 and "." in sys.argv[2]: # set a different extension
                    extension = sys.argv[2]
                filelist = os.listdir(sys.argv[1])
                dir=sys.argv[1]
            
        elif optioncount == 2:
            if "." in sys.argv[1]:
                infile = sys.argv[1]
                filelist.append(infile)
            else:
                dir = sys.argv[1]
                filelist = os.listdir(sys.argv[1])

        outfolder = ".\\"                          # Assume user has specified a local file
        if infile.rfind("\\") > 1:                 # unless a folder is in the name.
            outfolder = infile[:infile.rfind("\\")]


        for n in range(len(filelist)):
            if extension in filelist[n]:
                if len(filelist)>1:
                    infile = dir+"/"+filelist[n]

                if infile.find('.') > 0:
                    outfile = infile[:infile.find('.')]+'.sac'
                    seedfile = infile[:infile.find('.')]+'.mseed'
                    
                else:
                    outfile = infile +'.sac'
                    seedfile = infile + '.mseed'
                outfolder = ".\\"                          # Assume user has specified a local file
                if infile.rfind("\\") > 1:                 # unless a folder is in the name.
                    outfolder = infile[:infile.rfind("\\")]

                PNE       = load(infile)


#                        PNE[0] is the header where:
#                        PNE[0][0][1] = The comment descriptor of the data
#                        PNE[0][1][1] = The station Network code
#                        PNE[0][2][1] = The station name
#                        PNE[0][3][1] = Component axis(Z,N, or E)
#                        PNE[0][4][1] = Reference time from Wavetrack
#                        PNE[0][5][1] = manually calculated start time
#                        PNE[0][6][1] = Calibration factor
#                        PNE[0][7][1] = Time correction in seconds

#                        Process PNE[0] differently now, so it's not related to order, just in case
#                        something (like network) ends up blank or missing.

                metadata  = {'comment':'','network':'','station':'','component':'','reftime':'01_JAN_1970_00:00:00.000', \
                             'starttime':'01_JAN_1970_00:00:00.000','cf':1.0,'tc':0.0}
                for data in PNE[0]: # Load the dictionary with the stuff in the header
                    if data:
                        metadata[data[0].lower()] = data[1]                
                Comment   = metadata['comment'] 
                Network   = metadata['network'] 
                Stname    = metadata['station'][:7] 
                Component = metadata['component'][:3] 
                CF        = np.float32(metadata['cf']) # Amplification of the system (crest factor) (now unused)
                first_sampletime    = float(PNE[1][0][0])# Number of seconds ahead of the ref_time in which the first sample actually occurs
#                                                          In New digitizations, first_sampletime is expected to be zero.
                Ref_time  = time.strptime(metadata['reftime'][:-4],"%d_%b_%Y_%H:%M:%S") # Tuple. This is apparent time of 1st digitized sample from the seismogram.
                Ref_time_Frac_sec = float(metadata['reftime'][20:])
                Start_time   = time.strptime(metadata['starttime'][:-4],"%d_%b_%Y_%H:%M:%S") # Tuple. This is the actual GMT time where the file will be sliced.
                Start_time_Frac_sec  = float(metadata['starttime'][20:])
                TC        = float(metadata['tc'])               # time correction constant from the seismogram that corrects apparent time to actual GMT


                reftime = calendar.timegm(Ref_time)+Ref_time_Frac_sec+float(first_sampletime) # Add the offset from the first sample and TC to reftime.
                starttime = calendar.timegm(Start_time)+Start_time_Frac_sec+float(first_sampletime) # Add the offset from the first sample to starttime.
                St_time = time.gmtime(starttime+TC)   # This is a tuple representing the corrected start time in GMT
                DeltaT = starttime - reftime  # time (in seconds) representing the data to be discarded at the beginning of stream  
                Starting_sample = float(DeltaT + first_sampletime)


#
#                       Delta is the average sample period. It is calculated from the offset time of last sample 
#                       minus offset of first sample / total number of samples.

            
                Delta     = (float(PNE[1][len(PNE[1])-1][0])-first_sampletime)/(len(PNE[1])-1)

#                       Load Data array
#                       Samples in file are (no longer) multiplied by 10,000 to convert from
#                       measurements of centimeters to microns, then it's divided by
#                       the Amplification (conversion) factor, known as CF (This is done within the dataless SEED file now)


                PNEtime = []
                PNEamplitude = []
#
#                                    Starting_sample is the index time representing the requested start time
#                                    Any samples that precede this index time are NOT added to the stream.
#
                for data in PNE[1]:
                    if float(data[0]) >= (Starting_sample):  # Discard the unwanted samples from beginning of file 
                        PNEtime.append(float(data[0]))
                        PNEamplitude.append(float(data[1]))

#
#                Process file for click points, then interpolate into smoother waveform

                PNEamp   = np.arange(len(PNEamplitude),dtype=np.float32)
                for n in range(len(PNEamplitude)):      
                    PNEamp[n] = np.float32(PNEamplitude[n]) 
                PNEamp = polynomial(PNEamp, order = 1, plot = False) # polynomial linear trend remove the offset.
 
                PNE_time,PNE_amp,PNEclicktime,PNEclickamp,errtime = Pchip(PNEtime,PNEamp)

                b        = np.arange(len(PNE_amp),dtype=np.float32)
                for n in range(len(PNE_amp)):                        #   Load the array with time-history data
                    b[n] = np.float32(PNE_amp[n])


                print( "Sample rate = {0:2.3f} \nNumber of samples = {1} \nFirst_sampletime = {2} seconds, last_sampletime = {3:2.3f} seconds".format( \
                       1/Delta,len(b),first_sampletime,float(PNE[1][len(PNE[1])-1][0])))
                print(f'Starting sample for the time slice occurs at {Starting_sample} seconds.') 
#
#                                            Create the SAC stream
#

                t        = SACTrace(data = b)         
                                             # set the SAC header values
                t.scale  = 1.0               # Set the scale for use with DIMAS software
                t.delta  = Delta
                t.nzyear = St_time.tm_year
                t.nzjday = St_time.tm_yday
                t.nzhour = St_time.tm_hour
                t.nzmin  = St_time.tm_min
                t.nzsec  = St_time.tm_sec
                t.nzmsec = int(metadata['starttime'][21:])          # int((Frac_second)*1000)
                t.kstnm  = Stname[:7]
                t.kcmpnm = Component
                t.IDEP   = 2                 # 4 = units of velocity (in Volts)
                                             # Dependent variable choices: 
                                             # (1)unknown, 
                                             # (2)displacement(nm), 
                                             # (3)velocity(nm/sec), 
                                             # (4)velocity(volts), 
                                             # (5)nm/sec/sec
                t.kinst  = "Displace"        # Instrument type
                t.knetwk = Network #         # Network designator
                t.kuser0 = "CM"        # Centimeters. Place the system of units into the user text field 0

                filename = Stname[:7]+"."+Network+".."+Component+"."+time.strftime('%Y.%m.%d.%H.%M.%S', St_time)
                outfile = os.path.join(outfolder,filename+".sac")
                with open(outfile,'wb') as sacfile:
                    t.write(sacfile)
                print (" File successfully written: {0}".format(outfile))       
                sacfile.close()


#                stats = {'network': Network, 'station': Stname[:7], 'location': '',\
#                     'channel': Component, 'npts': len(b), 'sampling_rate': 1/Delta,\
#                     'mseed': {'dataquality': 'D'}} 
#                stats['starttime'] = UTCDateTime(mkgmtime(St_time)+Frac_sec)

#                st=Stream([Trace(data=b, header = stats)])
#                st[0].data = st[0].data.astype('int32')  # dtype it for use with steim2
#                st[0].write(seedfile,format="MSEED", encoding = 'STEIM2')

                st=read(outfile)
                seedfile = os.path.join(outfolder,filename+".mseed")
                st.write(seedfile,format="mseed")
                print (" File successfully written: {0}".format(seedfile))
                print ("File written to {}".format(outfile))

                plt.figure(figsize=(15,5))
                plt.title(f'Wavetrack vs. PCHIP plot for {filename}')
                plt.xlabel(f'Time in seconds, at {1/Delta:.2f} samples/sec')
                plt.ylabel(' Amplitude (in centimeters) ')
                plt.scatter(PNEtime, PNEamp, s=1, c='blue',label = "Original Wavetrack Output", alpha=0.35)
                plt.scatter(PNEclicktime,PNEclickamp, s=5, c='green', label = "Extracted click points", alpha = 1.0)
                plt.plot(PNE_time, b,c='red',label = "PCHIP interpolation")
                plt.bar(errtime, height = 5, width=.02, bottom=-2.5, align='center',color = "green",fill=False)
                plt.legend()
                plt.show()
        
    else:
        print ("Useage: PNE2SAC infile.txt (outfile.asc)")
        print ("Or, PNE2SAC target_directory target_extension(like .txt)")
        print ("No infile or directory specified.")
        print (len(sys.argv))

#
# Check and run the main function here:
#
if __name__ == '__main__':
  main()
