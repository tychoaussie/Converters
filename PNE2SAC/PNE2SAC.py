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
__version__ = "20210519"
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

import os, csv, sys, numpy as np #, time, calendar, string
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
# version 20210519: Fix the SAC timing that was missing the time correction constant.
# version 20210428: Revision 2.0 includes:
#                   Waveform analysis and optimization algorithm, support of clickpoint input files
#                   Revised plot of before/after, improved miniseed file generation
#                   Modify how time is handled, and move to UTCDateTime and convert Time_params to a DICT.
# version 20201205: Include SAC header field 'lpspol' to explicitly set polarity of the channel
# Version 20201119: Add the channel location to the export file so that alternate location codes can be added.
# Version 20201112: Include channel orientation fields in the SAC file to correctly convert to polar coordinate system in SAC
# Version 20200831: Accomodate dates earlier than 1970
#         20200813: Polarity set to 1.0 in this code
# Version 20200730: Add minor bug fix in def load() for when incoming text files contain blank lines with spaces inside.
# Version 20200829: Add a total to the number of discontinuities found in the wavetrack output that need
#                   correction, and place this information into the graph.
# Version 20200504: Fixed the start time fraction of second Frac_sec which was incorrect.
#                   Also, the graph now gets save by default. Polarity has been corrected, as well. 
#                   Future changes:  The text log will also be saved by default,
#                   and filename will be time coded with the time of miniseed generation.
#                   Make the SAC file format an OPTION, and only output it when directed with a '-SAC' switch. 
#                   Add Polarity awareness. Does a negative-going wavetrack measurement represent a positive
#                   or a negative-going waveform? In any event, make lower-going values reflect POSITIVE as default.
# Version 20200318: Add in the Mouseclick file output that provides the interpolated mouse clicks.
# Version 20200217: Minor changes to how the start time is generated (specifically the first sample offset)
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
# Each text file must contain a header with the following nine rows (comment and location are optional):
#
#    COMMENT Neva-2-2-Stolb-409_600dpi
#    NETWORK RY
#    STATION STB
#    COMPONENT HHZ
#    LOCATION (optional). Default is "blank" but other common values are 00,01
#    REFTIME   24_JUL_1987_02:03:00.000
#    STARTTIME 24_JUL_1987_02:05:00.000
#    CF 25000
#    TC -1.5
#    OPTIMIZE TRUE (optional). False is the default.
#        You can set it TRUE if you want optimization algorithm to make changes automatically.
#    THRESHOLD 8 (optional. 8 is the default for optimization. This can be set between 4 and around 12, 
#        where lower numbers are more aggressive with the optimization. 
#        Optimize moves clickpoints back & forth to minimize glitches in the
#        derivative above 8 Hz but if you are too aggressive it makes things worse. 
#
#    Headers are followed by the data fields, representing time (relative to REFTIME) and amplitude (in centimeters):
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
# version 20201205: Include SAC header field 'lpspol' to explicitly set polarity of the channel

#                                   The Defs for PNE2SAC:
#
def str2bool(v):                         # Return a boolean value if text is found that matches the list
  return v.lower() in ("yes", "y", "true", "t", "1")

# FIND_SPIKES: Find locations in the data where input is the list of data from the 1-pole 8hz filter stream, 
# corresponding PNE_time in seconds since first sample, and the threshold vale
# Output is the sample time of the spike in seconds since first sample, and it's corresponding value
#
def Find_spikes(data,PNE_time,threshold):
    absdata = np.abs(data)
    samplenumber = []
    value = threshold*np.mean(absdata)   # Trigger when value exceeds n times the mean
    for i,sample in enumerate(absdata):
        if sample > value:
            samplenumber.append(i)
#    print (len(samplenumber))

    sampletime = []
    amplitude = []
    for sample in samplenumber:
        sampletime.append(PNE_time[sample])
        amplitude.append(data[sample])
    return(sampletime,amplitude) # representing the points where threshold was exceeded
#
# Input is the spike's sample time in seconds since the first sample, and the list of clickpoint times.
# Output is the matched clickpoint in seconds since first sample, and the index number of said clickpoint in the list.
#

def Match(sampletime,clicktime): # delta is time in seconds that the sample must be nearest in order to be considered a match
    match = []
    match_index = []
    clickminimum = find_minimum_click_interval(clicktime)
    for sample in sampletime: # sampletime is the time where the sample exceeds threshold
        for click in clicktime: # match is the clicktime that preceeds to the exceeded threshold.
            if np.abs(click-sample) <= 2*clickminimum:
                if not click in match:
                    match.append(click)
                    match_index.append(clicktime.index(click))
                    #print("{0} {1:0.3f} {2:0.3f}".format(clicktime.index(click),click,clicktime[clicktime.index(click)]))
    print(f'Total number of matches: {len(match)}')
 #   print(len(sampletime))
    return(match,match_index) #returns the clickpoints that match up with samples that exceeded the trigger

#
#  Find the minimum interval within a clickfile series, that represents probably a single pixel of resolution.
#

def find_minimum_click_interval(clicktime):
    clickminimum = clicktime[len(clicktime)-1]-clicktime[0]
    for i in range(1,len(clicktime)):
        interval = clicktime[i]-clicktime[i-1]
        if interval < clickminimum:
            clickminimum = interval
    return(clickminimum)
#
# Input is the clickpoint list, along with the index number of the selected clickpoint in need of adjustment.
# Output is the refined clickpoint list.
#
def Adjust(match_index,clicktime,shiftlimit):
    # Adjust the clicktimes associated with the glitch to reduce the excessive velocities
    clickminimum = find_minimum_click_interval(clicktime)  # Estimated Time interval for one pixel 
    for index in match_index:
        if index > 0: # Beware of the penultimate sample in original
            if index < (len(clicktime)-1):
                original = clicktime[index]
                pre = clicktime[index-1]-clicktime[index]
                post = clicktime[index+1]-clicktime[index]
                # One thing: assume the point is off by only a few pixels.
                # Shift, but limit it to maximum pixel shift of 'n' pixels
                if shiftlimit < clickminimum: # Dont let the shiftlimit fall below the one-pixel threshold
                    shiftlimit = clickminimum 
                maxshift = shiftlimit # *clickminimum # max time in seconds allowed in a given clickpoint adjustment
                # Preserve polarity
                polarity = ((post+pre)/2.0)/np.abs((post+pre)/2.0)

                if np.abs((post+pre)/2.0) < maxshift:
                    shift = (post+pre)/2.0
                else:
                    shift = polarity*maxshift

                clicktime[index] = clicktime[index]+shift
                print(f'Adjusting clickpoint[{index}] at {original:0.3f} by {(clicktime[index]-original):0.3f} sec. to {clicktime[index]:0.3f} seconds.')
            else:
                print(f"Last clickpoint, although tagged for adjustment, has not been modified.")
    return(clicktime) # return a modified set of click points

#
# Optimize function checks for unrealistic spikes and optimizes clickpoint placement to minimize high-velocity
# spikes above the assumed seismometer channel operating bandwidth
# then provides a modified list of clickpoints along with diagnostic lists
# such as the clickpoints adjusted and related amplitudes,
# the testdata derivative, and the spikes that exceed the threshold

def Optimize(clicktime,clickamp,streamtime,streamdata,metadata,time_params):
    # This is a hard rule for now, but may change.
    #
    tr = Maketrace(streamdata,metadata,time_params) # makes a trace from streamdata, a time-history output from PCHIP
    #print(tr.stats)
    tr2 = tr.copy()
    tr2.differentiate(edge_order=2)
    tr2.filter("highpass",freq = metadata['breakpoint'], corners = 2) # single-pole filter above operating frequency to find big spikes
#    tr2.differentiate()

    testdata = [(sample) for sample in tr2.data] # convert it to an absolute value as we look for excursions

    spiketime,spikeamp = Find_spikes((testdata),streamtime,float(metadata['threshold'])) # test against a threshold
    # Match requires the spiketime match the clicktime within +- delta
    #
    # find minimum time between click points
    #
    clickminimum = find_minimum_click_interval(clicktime)
    print(f'minimum time between clicks is {clickminimum:0.3f} seconds.')

    Clicktime_to_adjust,Clicktime_to_adjust_index = Match(spiketime,clicktime) # find problem clickpoints that made the spikes

#    print(len(Clicktime_to_adjust),len(Clicktime_to_adjust_index))

    new_clicktime = Adjust(Clicktime_to_adjust_index,clicktime,metadata['shiftlimit']) # Make revised clickpoints

    C2amp = []                   # Make an adjusted Clickpoint_to_amplitude list
    for index in Clicktime_to_adjust_index:
        C2amp.append((clickamp[index]))   

    # Optimize will provide the following:
    # Revised list of clicktimes, and amplitudes,
    # the filtered trace for plotting,spiketime,spikeamp,
    # list of original clickpoint times, and the list of amplitudes associated with those original points.
    return(new_clicktime,clickamp,Clicktime_to_adjust,C2amp,testdata,spiketime,spikeamp)


#
# Import clickfile imports the previously generated clickfile and outputs the dictionary header, along with
# three lists representing the absolute time, relative time, and amplitude of the clickpoint originally
# extracted from a wavetrack text file.
# Output is header with necessary metadata
# absolute time is since first sample

def Import_clickfile(clickfile):
    with open(clickfile,'r') as fin:
        list = csv.reader(fin)
        clicktime = []
        clickamp = []
        metadata = {'comment':'','station':'','network':'','component':'','reftime':'01_JAN_1970_00:00:00.000', \
        'starttime':'01_JAN_1970_00:00:00.000','cf':float(1.0),'tc':float(0.0), 'location':'', 'samplerate':float(100.0),\
        'nsamples':int(0),'polarity':1.0,'threshold':16.0,'breakpoint':16.0,'optimize':'False','shiftlimit':0.085}
        for row in list:
            if len(row) > 0:          # skip blank lines
                r = row[0].split()    # First seven rows are header and assigned to PNE[0]
                if len(row):
                    if  row[0].lower() in metadata.keys():
                        metadata[row[0].lower()] = row[1]   
                        print(f' Dict item {row[0]} set to {row[1]}')
                    else:                 # Stack represents data and is assigned to PNE[1]
                        try:
                            clicktime.append(float(row[1]))
                            clickamp.append(float(row[2]))
                        except:
                            print(f'Ignoring import of row{row}')
        # Type the values for cf, tc, samplerate, and nsamples
    metadata['breakpoint'] = float(metadata['breakpoint'])
    metadata['shiftlimit'] = float(metadata['shiftlimit'])       
    metadata['cf'] = float(metadata['cf'])
    metadata['tc'] = float(metadata['tc'])
    metadata['samplerate'] = float(metadata['samplerate'])
    metadata['nsamples'] = int(((clicktime[len(clicktime)-1]-clicktime[0])*metadata['samplerate'])+1)
    metadata['optimize'] = str2bool(metadata['optimize'])

    first_sampletime    = float(clicktime[0])# Number of seconds ahead of the ref_time in which the first sample actually occurs
    Time_params = Get_time_params(metadata,first_sampletime)

    #
    #                       Delta is the average sample period. It is calculated from the offset time of last sample 
    #                       minus offset of first sample / total number of samples.
    # 

    metadata['delta']     = (float(clicktime[len(clicktime)-1])-first_sampletime)/(metadata['nsamples'])
    print(f"delta has been set to {metadata['delta']} seconds, which is {1/metadata['delta']} Hz")
    errtime = False
    return(clicktime,clickamp,errtime,metadata,Time_params)


def Export_clickfile(metadata,clicktime,clickamp):
#
#               Write the mouse click file.
#
    clickfile = os.path.join(metadata['outfolder'],(metadata['filename']+"_clickpoints.csv"))
    with open(clickfile, mode='w', newline='\n') as mouseclick:
        mouseclick_writer = csv.writer(mouseclick, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        mouseclick_writer.writerow(['COMMENT',metadata['comment']])
        mouseclick_writer.writerow(['STATION',metadata['station']])
        mouseclick_writer.writerow(['NETWORK',metadata['network']])
        mouseclick_writer.writerow(['COMPONENT',metadata['component']])
        mouseclick_writer.writerow(['CF',metadata['cf']])
        mouseclick_writer.writerow(['TC',metadata['tc']])
        mouseclick_writer.writerow(['REFTIME',metadata['reftime']])
        mouseclick_writer.writerow(['STARTTIME',metadata['starttime']])
        mouseclick_writer.writerow(['OPTIMIZE',metadata['optimize']])
        mouseclick_writer.writerow(['THRESHOLD',metadata['threshold']])
        mouseclick_writer.writerow(['SHIFTLIMIT',metadata['shiftlimit']])
        mouseclick_writer.writerow(['BREAKPOINT',metadata['breakpoint']])
        mouseclick_writer.writerow([])
        mouseclick_writer.writerow(['Clickpoint','Relative_time(seconds)','Amplitude(centimeters)'])
        for i, mouseclicks in enumerate(clicktime,0):
            mouseclick_writer.writerow([i,f'{(clicktime[i]-clicktime[0]):003.003f}',clickamp[i]])
    print ("Mouseclick file successfully written: {0}".format(clickfile))

#
# Load a wavetrack output file with header.
# header contains the necessary metadata
# stack contains two element list, of time in seconds since first sample, amplitude in centimeters.
#    
def Import_Wavetrack(infile):
    with open(infile,'r') as fin:
        list = csv.reader(fin)
        PNE = []
        metadata = {'comment':'','station':'','network':'','component':'','reftime':'01_JAN_1970_00:00:00.000', \
        'starttime':'01_JAN_1970_00:00:00.000','cf':float(1.0),'tc':float(0.0), 'location':'', 'samplerate':float(100.0),\
        'nsamples':int(0),'polarity':1.0,'threshold':16.0,'optimize':'False','breakpoint':16.0,'filename':'none',\
        'shiftlimit':0.085,'outfolder':'C:\\'}
        for row in list:
            if len(row) > 0:          # skip blank lines
                r = row[0].split()    # First seven rows are header and assigned to PNE[0]
                if len(r) > 1:
                    if  r[0].lower() in metadata.keys():
                        metadata[r[0].lower()] = r[1]    # rowcnt 
                        print(f' Dict item {r[0]} set to {r[1]}')
                    else: # Stack represents data and is assigned to PNE[1]
                        try:
                            PNE.append(r)
                        except:
                            print(f'Ignoring row {r}')
#
#                           Make sure the numerical constants in metadata are interpreted as numbers
#
    metadata['nsamples']   = int(metadata['nsamples'])
    metadata['polarity']   = float(metadata['polarity'])
    metadata['threshold']  = float(metadata['threshold'])
    metadata['breakpoint'] = float(metadata['breakpoint'])
    metadata['samplerate'] = float(metadata['samplerate'])
    metadata['cf']         = float(metadata['cf'])
    metadata['tc']         = float(metadata['tc'])
    metadata['shiftlimit'] = float(metadata['shiftlimit']) # shiftlimit is listed in seconds.
    metadata['optimize'] = str2bool(metadata['optimize'])   # Boolean value based in interpretation of text field
    #
    #                       Build the time parameters dictionary
    #

    first_sampletime    = float(PNE[0][0])# Number of seconds ahead of the ref_time in which the first sample actually occurs 
    Time_params = Get_time_params(metadata,first_sampletime)

    #                       In New digitizations, first_sampletime is expected to be zero.
#    print(f"\nStarttime is calculated as: {Time_params['starttime']}")
    #
    #                       Delta is the average sample period. It is calculated from the offset time of last sample 
    #                       minus offset of first sample / total number of samples.

    Delta     = (float(PNE[len(PNE)-1][0])-first_sampletime)/(len(PNE)-1)
    metadata['delta'] = Delta                # This is the calculated sample interval, which should coincide with the stated sample rate
    metadata['nsamples'] = len(PNE)          # remember to adjust this value when trimming the file to size
    return(metadata,PNE,Time_params)



#
#
# Get Time Parameters: Parse out the start time and reference time, generate UTCDateTime instances
# Times listed for Starttime and Reftime are pre-time correction
# And compute the difference between start time and reference time, which is starting_sample in seconds
# beyond the beginning of the raw digitization. This is where we'll make the cut when generating the output.
# Add time correction TC AFTER making the cut, which is done later in the program.  
#
def Get_time_params(metadata,first_sampletime):
    TC        = float(metadata['tc'])   # time correction constant that corrects apparent time to actual GMT
    starttime = metadata['starttime']
    reftime = metadata['reftime']
    month = {'JAN':'01','FEB':'02','MAR':'03','APR':'04','MAY':'05','JUN':'06','JUL':'07','AUG':'08','SEP':'09', \
             'OCT':'10','NOV':'11','DEC':'12'}
    sttime = "%s-%s-%sT%s:%s:%s" % (starttime[7:11],month[starttime[3:6]],starttime[0:2],starttime[12:14],starttime[15:17],starttime[18:24])
    St_time = UTCDateTime(sttime)
    rftime = "%s-%s-%sT%s:%s:%s" % (reftime[7:11],month[reftime[3:6]],reftime[0:2],reftime[12:14],reftime[15:17],reftime[18:24])
    Rf_time = UTCDateTime(rftime)    
#    St_time = St_time+TC   # This is a tuple representing the corrected start time in GMT
    DeltaT = St_time-Rf_time  # time (in seconds) representing the data to be discarded at the beginning of stream  
    Starting_sample = float(DeltaT) # + first_sampletime)
    Time_params = {}
    Time_params["St_time"] = St_time # [0]
    Time_params["Rf_time"] = Rf_time
    Time_params['Starting_sample'] = Starting_sample # [1]
    Time_params['starttime'] = starttime # [2]
    Time_params['reftime'] = reftime # [3]
    return(Time_params)

#                       Load Data array from the Wavetrack output file
#                       Samples in file are (no longer) multiplied by 10,000 to convert from
#                       measurements of centimeters to microns, then it's divided by
#                       the Amplification (conversion) factor, known as CF (This is done within the dataless SEED file now)
#                       Use PCHIP to provide a list of click times, and interpolated waveform from those click times.

def PNE2Clicktime(metadata,PNE,Time_params):

    PNEtime = []
    PNEamplitude = []
    #
    #                                    Starting_sample is the index time representing the requested start time
    #                                    Any samples that precede this index time are NOT added to the stream.
    #
    for data in PNE:
        if float(data[0]) >= (Time_params['Starting_sample']):  # Starting_sample:  # Discard the unwanted samples from beginning of file 
            PNEtime.append(float(data[0])) # This is the relative time
            PNEamplitude.append(2-float(data[1])*metadata['polarity']) # This is amplitude in centimeters

    #
    #                Process file for click points, then interpolate into smoother waveform

    metadata['nsamples'] = len(PNEtime)    # Re-adjust the metadata to account for the trimming of the excess samples
    PNEamp   = np.arange(len(PNEamplitude),dtype=np.float32)
    for n in range(len(PNEamplitude)):      
        PNEamp[n] = np.float32(PNEamplitude[n]) 
    PNEamp = polynomial(PNEamp, order = 1, plot = False) # polynomial linear trend remove the offset.

    # Get the clickpoints

    PNEclicktime,PNEclickamp,errtime = WT_to_clickpoint(PNEtime,PNEamp)

    PNE_time,PNE_amp = Pchip(PNEclicktime,PNEclickamp,len(PNEtime))
    #
    #               Prepare a list of UTC corrected times that correspond to the click points.
    # PNEclicktime is a list of clickpoints in UTCDateTime
    PNE_clicktime = [Time_params['St_time']+(sample - Time_params['Starting_sample'])+float(metadata['tc']) for sample in PNEclicktime] 
    Clicktime = []
    for click in PNE_clicktime:
        Clicktime.append(str(click))
    return(PNE_time,PNE_amp,PNEclicktime,PNEclickamp,Clicktime,errtime)

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

def WT_to_clickpoint(time,amp): # input the time and amplitude series. method describes the desired interpolation method
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
    return(clicktime,clickamp,errtime)

def Pchip(clicktime,clickamp,nsamples): # Bring in clickpoints and output a PCHIPped time series
    # At this point, clicktime represents the time of each click point and clickamp represents amplitude of each mouse click
   # print(f'length of clicktime = {len(clicktime)} Len clickamp = {len(clickamp)} Len slopes = {nsamples}')
    x_data = np.array(clicktime)
    y_data = np.array(clickamp)

    Pchip_time = np.linspace(min(x_data), max(x_data), nsamples) # From the oldest time to the youngest time, space out 'n' points, as specified in len(slopes)
    bi = interpolate.PchipInterpolator(x_data, y_data)
    Pchip_amplitude = bi(Pchip_time)
    return(Pchip_time,Pchip_amplitude)

def Create_sacstream(b,metadata,time_params,outfolder):
    #                                          Create the SAC stream
    St_time = time_params['St_time']  + float(metadata['tc'])          # Shortcut for filename below
    #
    Delta     = (1/float(metadata['samplerate']))
    Network   = metadata['network']
    Component = metadata['component']
    Location  = metadata['location']
    Stname    = metadata['station']
    t         = SACTrace(data = b)         
                                 # set the SAC header values
    t.scale  = 1.0               # Set the scale for use with DIMAS software
    t.delta  = Delta
    t.nzyear = St_time.year
    t.nzjday = St_time.julday
    t.nzhour = St_time.hour
    t.nzmin  = St_time.minute
    t.nzsec  = St_time.second
    t.nzmsec = St_time.microsecond/1000          # int((Frac_second)*1000)
    t.kstnm  = metadata['station'][:7]
    t.kcmpnm = metadata['component']
    t.khole  = metadata['location']
    #                t.lpspol = True # True for left-hand-rule?
    Orientation = {'z':[0,0],'n':[0,90],'e':[90,90],'1':[0,90],'2':[90,90]} # Orientation of the channel must be included for SAC users
    t.cmpaz  = Orientation[str.lower(Component)[2:]][0]
    t.cmpinc = Orientation[str.lower(Component)[2:]][1]
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
    #

    ftime = (f"{St_time.year}.{St_time.month:02.0f}.{St_time.day:02.0f}.{St_time.hour:02.0f}.{St_time.minute:02.0f}.{St_time.second:02.0f}")
    filename = metadata['network']+"."+metadata['station']+"."+metadata['location']+"."+metadata['component']+"."+ftime

    outfile = os.path.join(outfolder,filename+".sac")

    #
    #               Write the SAC file
    #
    with open(outfile,'wb') as sacfile:
        t.write(sacfile)
    print (" File successfully written: {0} with a start time of {1}".format(outfile,St_time))       
    sacfile.close()
    return(t) #  This is the trace that was just written to disk.

def Make_obspy_trace(tr,metadata): # This converts the sac trace tr into an obspy trace then returns an obspy trace
    tr2 = tr.to_obspy_trace()
    #tr.stats.comment    = metadata['comment'] 
    tr2.stats.network   = metadata['network'] 
    tr2.stats.station   = metadata['station'][:7] 
    tr2.stats.channel   = metadata['component'][:3]
    tr2.stats.location  = metadata['location']
    tr2.stats.starttime = starttime
    tr2.stats.delta = 1/float(metadata['samplerate'])
    tr2.stats.npts= len(tr2.data)
    tr2.stats.calib = metadata['cf']
    tr2.stats
    return(tr2)
#
# Maketrace generates an Obspy trace suitable for writing out as a miniseed file
# as well as useful for plotting, and applying Obspy signal analysis tools.
#

def Maketrace(streamdata,metadata,time_params):
    tr = Trace(data=streamdata)
    tr.data = tr.data.astype('float32')
    tr.stats['delta']=metadata['delta']
    tr.stats['network'] = metadata['network']
    tr.stats['station'] = metadata['station'][:7]
    tr.stats['location'] = metadata['location']
    tr.stats['channel'] = metadata['component'][:3]
    tr.stats['starttime'] = time_params['St_time']+float(metadata['tc']) # Add the time correction here
    return(tr)



def PNEplot(PNE,clicktime,clickamp,PNE_time,streamdata,errtime,metadata,time_params):
    PNEtime = []
    PNEamplitude = []
    if len(PNE):
        # Do this stuff for plotting of the wavetrack data
        # clicktime,clickamp are from the clickpoint output
        # PNE_time,streamdata are from the PCHIP algorithm
        # errtime (if any) are clickpoints with discontinuities
        # Wavetrack channel is PNE[0],PNE[1]
        #
        for data in PNE:
            if float(data[0]) >= (time_params[5]):  # Starting_sample:  # Discard the unwanted samples from beginning of file 
                PNEtime.append(float(data[0])) # This is the relative time
                PNEamplitude.append(2-float(data[1])*metadata['polarity']) # This is amplitude in centimeters
        PNEamp   = np.arange(len(PNEamplitude),dtype=np.float32)
        for n in range(len(PNEamplitude)):      
            PNEamp[n] = np.float32(PNEamplitude[n]) 
        PNEamp = polynomial(PNEamp, order = 1, plot = False)       
    if(len(PNEtime)):
        plt.plot(PNEtime,PNEamp,color='black')
    if(len(PNE_time)):
        plt.plot(PNE_time,streamdata,color='green') # Show the PCHIP output
    plt.scatter(clicktime,clickamp,color='red',s=4) # Show the clickpoints
    plt.show()

# plot_optimized makes a plot with all the information, from non-optimized/optimized comparisons
# to clickpoints, to some stats about the digitization.
# It also calculates some statistics about the streams.


def plot_optimized(metadata,time_params,outfolder,tr1,tr2,errtime,streamtime,clicktime,Clicktime_to_adjust,spiketime,clickamp,spikeamp,testdata,newtestdata,C2amp):
    St_time = time_params['St_time']
    tr4 = tr2.copy()
    tr4.differentiate(edge_order=2)
    tr4.filter("highpass",freq = metadata['breakpoint'], corners = 2) # single-pole filter above operating frequency to find big spikes
#    tr4.differentiate()
    newtestdata = [sample for sample in tr4.data] # Newly corrected stream

    # Report whether or not the file being written has been modified with the optimization algorithm
    if metadata['optimize']:
        waveselect = "optimized"
    else:
        waveselect = "original"
    if errtime:
        errors = "Number of discontinuities in this file = "+str(len(errtime))
    else:
        errors = "Sourced from click file."
    adjusteds = "Number of clickpoints that exceed threshold = "+str(len(Clicktime_to_adjust))
    optimized = (f"Optimization algorithm: {metadata['threshold']} x rms trigger w/ {metadata['breakpoint']} Hz highpass filter.")
    shiftlimit = (f"Shift limit = {metadata['shiftlimit']}")
    outputselect = (f"Final output file represents the {waveselect} waveform.")
    fig, axs = plt.subplots(3)
    fig.suptitle(f"{optimized}\n{shiftlimit}\n{errors}\n{adjusteds}\n{outputselect}")
    fig.text(0.06, 0.5, 'Amplitude (centimeters)', ha='center', va='center', rotation='vertical')
    fig.set_size_inches(18.0,9.0)

    xmin = 0
    xmax = clicktime[len(clicktime)-1]
    #print(len(clicktime),len(Clicktime_to_adjust))
    if len(Clicktime_to_adjust) > 1:
            if Clicktime_to_adjust[0] > 10:
                xmin = Clicktime_to_adjust[0]-10
            if  (clicktime[len(clicktime)-1] - Clicktime_to_adjust[len(Clicktime_to_adjust)-1])  > 10: 
                xmax = Clicktime_to_adjust[len(Clicktime_to_adjust)-1]+10
    else:
        xmin = clicktime[0]
        xmax = clicktime[len(clicktime)-1]
#    xmin = 200.75
#    xmax = 241.0
    axs[0].set_xlim([xmin,xmax])
    axs[1].set_xlim([xmin,xmax])
    axs[2].set_xlim([xmin,xmax])
#    axs[2].set_ylim([-0.6,0.15])
#    axs[0].xlabel("original vs optimized")
    axs[0].annotate('Original(blue) vs optimized waveform(red)',xy=(xmin+(xmax-xmin)/2.0, np.min(tr2.data)))
    axs[1].annotate('optimized waveform (red)',xy=(xmin+(xmax-xmin)/2.0, np.min(tr2.data)))
    axs[2].annotate('Trigger channel before(blue) vs. after(red)',xy=(xmin+(xmax-xmin)/2.0, np.min(testdata)))
    axs[0].bar(errtime, height = np.max(tr1.data), width=.02, bottom=np.min(tr1.data), align='center',color = "green",fill=False)
    axs[0].plot(streamtime,tr1.data,label = "Original trace",color='blue',linewidth = 2) # original trace in blue
    axs[0].scatter(clicktime,clickamp,color = 'blue',s=8) # original clickpoints
    axs[0].plot(streamtime,tr2.data,color='red',label = "Modified trace",linewidth = 1) # optimized trace overlay
    axs[1].plot(streamtime,tr2.data,color='red') # optimized trace n second plot
    if C2amp:
        axs[0].scatter(Clicktime_to_adjust,C2amp,color = 'red',s=24) # adjusted clickpoints
    if spiketime:
        axs[2].scatter(spiketime,spikeamp,color = 'red',s=8) # Where threshold is exceeded
    if testdata:
        axs[2].plot(streamtime,testdata,color='blue', linewidth = 2) # trigger channel
    axs[2].plot(streamtime,newtestdata,color='red',linewidth = 1)
#    plt.text(xmin,np.min(tr1.data),errors)
    ftime = (f"{St_time.year}.{St_time.month:02.0f}.{St_time.day:02.0f}.{St_time.hour:02.0f}.{St_time.minute:02.0f}.{St_time.second:02.0f}")
    filename = metadata['network']+"."+metadata['station']+"."+metadata['location']+"."+metadata['component']+"."+ftime
    filenametitle = (f"Digitized channel name: {os.path.join(outfolder,filename)}.mseed\n{metadata['comment']}")
    fig.text(0.5, .0625, filenametitle, ha='center', va='center' , rotation='horizontal')
    plt.savefig(os.path.join(outfolder,filename+".png"))
    plt.show()

def main():
#                               Parse the command line switches
    optioncount = len(sys.argv)
    Folder = False
    SAC = True      # Always make sac files in this version.
    MSEED = True
    QCcheck = True
    Sps = 100
    outputfile_defined = False
    filelist = []
    Wavetrack = True
    Clickfile = False
#    metadata = {'comment':'','station':'','network':'','component':'','reftime':'01_JAN_1970_00:00:00.000', \
#        'starttime':'01_JAN_1970_00:00:00.000','cf':float(1.0),'tc':float(0.0), 'location':'', 'samplerate':float(100.0),\
#        'nsamples':int(0),'polarity':1.0,'threshold':8.0,'optimize':'False','shiftlimit':1.0}

    dir=""
    extension = '.txt'

    if optioncount > 1:
        for i,args in enumerate(sys.argv): # Check to see if alternate sample rate is specified
            if '-sps' in args.lower():
                Sps = float(args[4:])
                print(f'Sample rate requested = {Sps} samples per second')
            elif '-sac' in args.lower():
                SAC = True
            elif '-clickfile' in args.lower():
                Clickfile = True
                Wavetrack = False
                extension = ".csv"
            elif '-folder' in args.lower():
                try:
                    Folder = True # This means the first instance following folder represents the folder name.
                    dir = sys.argv[i+1]
                    filelist = os.listdir(sys.argv[i+1])
                except:
                    print("Warning! No folder specified within arguments. Check syntax.")
            elif '.' in args.lower() and i!=0: # A file name has been specified.
                infile = args
                filelist.append(infile)
            elif '-display' in args.lower():
                SAC = False
                MSEED = False
#        print(f"Clickfile = {Clickfile}, Wavetrack = {Wavetrack}, filename = {infile}")

        outfolder = ".\\"                          # Assume user has specified a local file
        if infile.rfind("\\") > 1:                 # unless a folder is in the name.
            outfolder = infile[:infile.rfind("\\")]

        for n in range(len(filelist)):
            if extension in filelist[n]:
                if len(filelist)>1:
                    infile = dir+"/"+filelist[n]

                #if infile.find('.') > 0:
                #    outfile = infile[:infile.find('.')]+'.sac'
                #    seedfile = infile[:infile.find('.')]+'.mseed'
                    
                #else:
                #    outfile = infile +'.sac'
                #    seedfile = infile + '.mseed'
                outfolder = ".\\"                          # Assume user has specified a local file
                if infile.rfind("\\") > 1:                 # unless a folder is in the name.
                    outfolder = infile[:infile.rfind("\\")]

#
#                 Import the wavetrack output file or the clickfile
#
                if Wavetrack:
                    metadata,PNE,time_params = Import_Wavetrack(infile)
                    streamtime,streamamp,clicktime,clickamp,GMTClicktime,errtime     = PNE2Clicktime(metadata,PNE,time_params)
                    streamdata        = np.arange(len(streamamp),dtype=np.float32)
                    for n in range(len(streamamp)): #   Load the array with time-history data
                        streamdata[n] = np.float32(streamamp[n])
                    St_time = time_params['St_time'] + float(metadata['tc'])
                    print(f"Importing Wavetrack file {infile}:")
                if Clickfile:
                    #print("Trying to import clickfile")
                    PNE = []
                    clicktime,clickamp,errtime,metadata,time_params = Import_clickfile(infile)
                    streamtime,streamamp     = Pchip(clicktime,clickamp,metadata['nsamples'])
                    streamdata        = np.arange(len(streamamp),dtype=np.float32) # b = time history intended for placing within the miniseed stream
                    for n in range(len(streamamp)):                                #   Load the array with time-history data
                        streamdata[n] = np.float32(streamamp[n])
                    St_time = time_params['St_time']+ float(metadata['tc'])
                    print(f"Importing Clickfile {infile}:")
                #
                #                             Print status lines 
                #
               # print(metadata)
                print( "\nSample rate = {0:2.3f} \nNumber of samples = {1}".format( \
                      1/float(metadata['delta']),len(streamdata)))
                if len(PNE):
                    print( "First_sampletime = {0} seconds, last_sampletime = {1:2.3f} seconds.".format( \
                           float(PNE[0][0]),float(PNE[len(PNE)-1][0])))
                    print(f"Starting sample for the time slice occurs at {float(time_params['St_time']-time_params['Rf_time'])} seconds.") 
#                print(f"St_time = {St_time.hour:02.0f}:{St_time.minute:02.0f}:{St_time.second:02.0f}.{int(St_time.microsecond/1000)}")
                print(f"Time correction = {metadata['tc']} seconds.")
                print(f"starttime = {time_params['starttime']}")
                if errtime:
                    print(f"There are {len(errtime)} discontinuities in this file.\n")
                    if len(errtime):
                        print("Repair discontinuities in Wavetrack output file, then re-run PNE2SAC.")
                #
                #              Calculate optimized trace
                #
#                if metadata['optimize']:
                new_clicktime,clickamp,Clicktime_to_adjust,C2amp,testdata,spiketime,spikeamp = \
                Optimize(clicktime,clickamp,streamtime,streamdata,metadata,time_params)

#                else:
#                    new_clicktime = clicktime
#                    Clicktime_to_adjust = [clicktime[0],clicktime[len(clicktime)-1]]
#                    C2amp = False
#                    spiketime = False
#                    testdata = False
#                    spikeamp = False
    
                New_streamtime,New_streamamp     = Pchip(new_clicktime,clickamp,metadata['nsamples']) # Make new PCHIP interpolation
                tr1 = Maketrace(streamdata,metadata,time_params) # makes a trace from old streamdata, a time-history output from PCHIP
                tr2 = Maketrace(New_streamamp,metadata,time_params) # generate a new trace from revised clickpoints.
                tr3 = tr1.copy()
                tr3.filter("highpass",freq = metadata['breakpoint'], corners = 1) # single-pole filter above operating frequency to find big spikes
                tr4 = tr2.copy()
                tr4.differentiate(edge_order=2)
                tr4.filter("highpass",freq = metadata['breakpoint'], corners = 2) # single-pole filter above operating frequency to find big spikes
#    tr4.differentiate()
                newtestdata = [sample for sample in tr4.data] # Newly corrected stream

                #
                # At this point, we have the original trace tr1, optimized trace tr2, and the trigger channel from the original at tr3.
                #
                #
                # Generate the file name for this dataset and place it in metadata to pass it to other defs for use in saving and plotting.
                #

                ftime = (f"{St_time.year}.{St_time.month:02.0f}.{St_time.day:02.0f}.{St_time.hour:02.0f}.{St_time.minute:02.0f}.{St_time.second:02.0f}")
                filename = metadata['network']+"."+metadata['station']+"."+metadata['location']+"."+metadata['component']+"."+ftime
                metadata['filename'] = filename
                metadata['outfolder'] = outfolder
                print(f"\nNumber of triggers in original channel = {len(spiketime)}")
                print(f"Number of clickpoints adjusted = {len(Clicktime_to_adjust)}")
                print(f"Maximum allowable shift of a clickpoint was defined as {metadata['shiftlimit']} seconds.")
                print(f"mean of original stream trigger channel = {np.mean(np.abs(testdata)):.2e} \nmean of modified stream trigger channel = {np.mean(np.abs(newtestdata)):.2e}")
                print(f"max spike amplitude of original stream = {np.max(np.abs(testdata)):.2e}\nmax spike amplitude of modified stream = {np.max(np.abs(newtestdata)):.2e}")
                print(f"%reduction in spike amplitude: {(1-np.max(np.abs(newtestdata))/np.max(np.abs(testdata)))*100.:02.1f} %\n")

                if QCcheck: # Assuming the waveforms passed quality check, save the files
                    #
                    # export the click file as either optimized or original
                    #
                    if SAC or MSEED:
                        if metadata['optimize']:
                            Export_clickfile(metadata,new_clicktime,clickamp)
                        else:
                            Export_clickfile(metadata,clicktime,clickamp)
                    #
                    # Write out the SAC file
                    #
                    if SAC:   # If user selected SAC files, write a SAC file
                        if metadata['optimize']:
                            print('Using the optimized trace to create the SAC file.')
                            Create_sacstream(New_streamamp,metadata,time_params,outfolder)
                        else:
                            print('Using the original trace without optimization to create the SAC file.')
                            Create_sacstream(streamdata,metadata,time_params,outfolder)
                    #
                    # Write out the miniseed file
                    #
                    if MSEED: # If user selected MSEED files, write a mseed file
                        outfile = os.path.join(outfolder,filename+".mseed")
                        if metadata['optimize']:
                            tr = tr2
                            print(f"Using the optimized trace to create miniseed file.")       
                        else:
                            tr = tr1 
                            print(f"Using the original trace to create miniseed file.")
                        tr.write(outfile,format = "MSEED")
                        print(f"File succesfully written: {outfile}")

                else: # Flagged QCcheck has failed so no output is generated.
                    print('Quality check failure; Check for discontinuities and re-edit the Wavetrack file.')
                    print('No SAC or miniseed files written to disk for this waveform.')        
                plot_optimized(metadata,time_params,outfolder,tr1,tr2,errtime,streamtime,clicktime,Clicktime_to_adjust,spiketime,clickamp,spikeamp,testdata,newtestdata,C2amp)
        
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