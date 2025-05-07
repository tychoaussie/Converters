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
__version__ = "20250408"
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

import argparse, time
import os, csv, sys, math, numpy as np
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
        
# Version 20250408: Correct SAC file outout to correctly calculate and apply the sample rate rather than use the default 100 sps, which is not always the case.
# Version 20250218: discontinuity error threshold in line_intercept() raised from .002 seconds to .010 seconds, corresponding to one sample at 100 samples/second.
#                 : Also change the threshold of WT_to_Clickpoint from 2 to 16 and save it in variable called dither_ratio to accomodate SKD instruments with long periods and low slope changes
# Version 20250113: Continued fillgaps() development. Address boundary conditions at first/last clickpoint of file if it lands on the minute mark. 
#                   Set velocity to zero on 1st and last clickpoint by setting it to the mean of the gap.
# version 20250110: output files are no longer created when there are discontinuities in the wavetrack export file. 
#                   PNE2SAC version number is output when code is started for tracking purposes.
#                   A new command line parsing method (argparse) has been added for more standardized program execution
#                   New options added to code including:
#                   -logfile for outputting the text to a log file without needing a pipe.
#                   -New fillgaps() improvements
#                   - Cleanup and enhanced documentation within the code
# version 20250107: Continued improvement of fillgaps plus integration of controls for automatic operation of code with python scripts from Jupyter notebook
#                   such as the -Nograph option. Use with Jupyer Notebook: "20250102_Parse_KAZAKH_PNE_folder_for_fillgaps_rebuild" for batch converting whole libraries
#                   from microsoft TEAMS containing digitalization libraries
# version 20240820: Add the optional fillgaps() to code for correcting timing marks
# version 20210810: Add alternate interpolation methods into the code for research purposes.
#                   New interpolation methods can be called by including additional
#                   header value in the input file. If left out, the default is PCHIP.
#                   Put in the header called 'interpolation' such as:
#                   interpolation pchip ( this is the default)
#                   interpolation cubicspline
#                   interpolation cubichermitespline
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
# Each text file must contain a header with the following rows (some are optional):
#
#    COMMENT Neva-2-2-Stolb-409_600dpi (optional)
#    NETWORK RY
#    STATION STB
#    COMPONENT EHZ
#    LOCATION (optional). Default is "blank" but other common values are 00,01
#    POLARITY (optional). Default is 1.0. Positive is left-going swings in the signal on wavetrack. Set to -1.0 if you suspect the positive-going wave faces to the right on the wavetrack image.
#    REFTIME   24_JUL_1987_02:03:00.000
#    STARTTIME 24_JUL_1987_02:05:00.000
#    CF 25000
#    TC -1.5
#    FILLGAPS TRUE (optional). False is the default. Set true if you want to fill the timing gaps with synthetic waveform. 
#    OPTIMIZE TRUE (optional). False is the default.
#        You can set it TRUE if you want optimization algorithm to make changes automatically.
#    THRESHOLD 25 (optional) (25 is the default. Most practical values range from 25 through 50)
#        where lower numbers are more aggressive with the optimization. 
#        Optimize moves clickpoints back & forth to minimize glitches in the
#        derivative above 16 Hz but if you are too aggressive it makes things worse.
#    BREAKPOINT 16 (optional). highpass filter corner frequency for derivative. 16 Hz is the default and seems to work for everything.
#    SHIFTLIMIT 0.085 (optional). This sets how many seconds the click point is allowed to move during optimization.
#        Shiftlimit values typically range between 0.085 and 0.250 for optimizing waveforms. 
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
#    and that the component match standard naming convention for IRIS/SEED channel names (i.e. EHZ,HLN,EHE, etc.)
#    Typically, a short-period sensor such as an SKM will use channel name EHZ or ELZ depending on gain
#    Use ELZ for gains less than about 1800. Use EHZ for gains higher than that.
#    Use the location code designator to separate out the mid-gain, normal gain, and high gain channels. 
#    Make it blank for normal gain short period channels, like channels from around 20000 to 50000
#    Set location code to 01 for mid-gain channels from 1800 through about 15000
#    Set location code to 02 for the high gain channels from 50000 to 100000+
#    This will depend on what you have for the station for that particular epoch. Some stations used multiple gains, some not.
#
#    Typically a longer period sensor such as an SKD (period over 10 seconds) will use HLZ,HLN,HLE
#


#                                   The Defs for PNE2SAC:
#

#
# New argument parser (08JAN2025) for program execution that looks more standard than my old hand-written code
#

parser=argparse.ArgumentParser(
    description=
    '''PNE2SAC is a utility for converting text files representing seismic waveforms.\n
       These text files are exported from Wavetrack. The text file is then modified by
       adding a header containing the metadata associated with the digitization.\n\n  
       Useage: Activate Python's ObsPy environment, then from the command prompt
       specify the target file.\n\n 
       Syntax: python PNE2SAC.py infile_name
       where, infile_name represents the output file from Wavetrack.\n\n
       Typical useage:\n
       <ObsPy> C:\\Python27\\scripts> python PNE2SAC.py C:/digitizations/station/channel/infile.txt\n\n
       Other options include:\n
       -display    : which opens and displays, but does not output any miniseed or SAC files\n
       -clickfile  : which opens a clickfile output from a previous PNE2SAC project rather than a wavetrack file.\n\n
       Consult the source code for more documentation on necessary header fields for the wavetrack output file and revision history.''',
    epilog="""Daniel Burk, Michigan State University mailto:burkdani@msu.edu""")
parser.add_argument('filename', help='PNE2SAC requires at minimum, the filename of the wavetrack output file with header added.')
parser.add_argument('-clickfile', action = 'store_true', help='Use a previously generated clickfile as the input file.')
parser.add_argument('-folder', action = 'store_true', help='input folder for conversion of multiple wavetrack output files at a time from the same folder')
parser.add_argument('-nograph',  action = 'store_true', help='Suppress graphical output to terminal and route text to log file.')
parser.add_argument('-logfile',  action = 'store_true', help='Route screen text to a log file in the destination folder.')
parser.add_argument('-optimize', action = 'store_true', help='Fine-tune clickpoint timing to minimize glitches in the derivative of the resulting signal.')
parser.add_argument('-fillgaps', action = 'store_true', help='Find 1-second timing gaps in data and fill them with synthetic signal to smooth for signal processing')
parser.add_argument('-display', action = 'store_true', help='Suppress file generation and only create the display graphic for viewing.')
parser.add_argument('-sps', type=float,default=100.0, help='Set an optional sample rate other than the default 100 samples/second on miniseed/sac output files')



def str2bool(v):                         # Return a boolean value if text is found that matches the list
  return v.lower() in ("yes", "y", "true", "t", "1")

# FIND_SPIKES: (used by Optimize) Find locations in the data where input is the list of data from the 1-pole 8hz filter stream, 
# corresponding PNE_time in seconds since first sample, and the threshold vale
# Output is the sample time of the spike in seconds since first sample, and it's corresponding value
# Discount clickpoints that are too small in terms of amplitude deviation to make a difference.
# Do this by comparing the difference in amplitude to the max/min of the 10 previous and 10 proceeding clickpoints
# If the difference in amplitude is less than 1/4 of the maxmin, don't tag it as a spike.
# 
def Find_spikes(data,PNE_time,threshold):
    absdata = np.abs(data) # Create the absolute value of the data stream, which is a list of amplitudes 
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
# def Match: (Used by Optimize)
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
    print(f'Total number of matches: {len(match)} out of {len(sampletime)}clickpoints exceed the trigger threshold specified.')
    return(match,match_index) #returns the clickpoints that match up with samples that exceeded the trigger

#
#  def find_minimum_clickpoints(Used by Optimize() and Adjust())
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
# Def Adjust, used by optimize():
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
                # Preserve polarity but watch out for when post+pre = 0
                polarity = 1
                if post+pre != 0: # Divide by the ABS of itself to capture polarity as being either 1.0 or -1.0 
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


# Used by fillgaps(): 
# generate_points : This is one of the frogDNA generating engines (10-16-2024).
# Note that the inflection points that bracket the timing mark do not represent true inflectin points, but instead represent where the timing gap
# begins and ends. Thus they must be dealt with as they introduce artifacts into the signal when removing instrumetn response.
# First synthesized clickpoint falls on the same line as incoming clickpoint to eliminate it's erroneous inflection point.
# Last synthesized clickpoint falls on the same line as the outgoing bracketing clickpoint as it also does not represent an inflection point
# Outputs from the def are a list of clickpoints, with synthesized time and amplitude for filling in the gap.


#  def generate_points3: Used by fillgaps(): 
#  This is the frogDNA generating engine.
#
# For scenarios in fillgaps, where the synthesized points are based on the starting and ending direction of the point amplitudes on initial velocities. 
# Four scenarios possible: (-V1 & -V2), (+V1 & +V2), (+V1 & -V2), (-V1 & +V2)
# Inputs to the def are:
    # generate_points inputs:
    # T1 = time of beginning point 
    # A1 = Amplitude of beginning point
    # V1 = Velocity at beginnign point
    # T2 = Time of endin point
    # A2 = Amplitude at ending point
    # V2 = Velocity at ending point
    # Numpoints = rate of points per second
    # tmed: median time gap
    # tmin: minimum time gap
    # Alow: lowest amplitude from adjacent clickpoints
    # Ahigh: highest amplitude from adjacent clickpoints
    # Outputs from the def are a list of clickpoints, with synthesized time and amplitude for filling in the gap.

        # generate_points3 = Build four clickpoints, with amplitudes that move from segments[2] to segments[4], but
        #                    with amplitudes of 75% peak from cpmax to cpmin.
        #                    Clockpoint[0] occurs at T1+(T2-T1/6) * 1 with amplitude: ((T2-T1)/6)* (i+1)
        #                    Clickpoint[1] occurs at T1+(T2-T1/6) * 2 with amplitude: (cpmax-cpmin)/2+(segments[3]-segments[2]/6) * 2
        #                    Clickpoint[2] occurs at T1+(T2-T1/6) * 3 with amplitude: -(cpmax-cpmin)/2+(segments[3]-segments[2]/6) * 3
        #                    Clickpoint[3] occurs at T1+(T2-T1/6) * 4 with amplitude: (cpmax-cpmin)/2+(segments[3]-segments[2]/6) * 4
        #                    Clickpoint[4] occurs at T1+(T2-T1/6) * 5 with amplitude: -(segments[3]-segments[2]/6)* 5
        # generate_points builds four clickpoints that are aware of the mean in the signal and add it to a pseudo-random signal.
        # Still some issues with conditions, where the initial & final amplitude going into the gap seem to have an effect.
        # Thus , point amplitudes need to be chosen to negate the apparent mean change over the timing gap such that it arrives
        # at zero displacement over the gap length.
        # Pass to the def the following : 
        # A1/V1/T1: Clickpoint & it's instantaneous velocity going into the gap
        # A2/V2/T2: Clickpoint & it's instantaneous velocity coming out of the gap
        # Numpoints
        # Tmed/Tmin/Alow/Ahigh/Amean = Stats about the clickpoints surrounding the gap such as variation in clickpoint density, and variation of the amplitudes.
        # One variable we might need is the mean on either side of the gap, for we want a "zero" total displacement across the gap (Amean).
        # It is expected however, that ultimate success cannot be achieved unless the code be made aware of the amplitude/phase curve of the channel and that
        # cannot be applied at the clickpoint level.
        # Original alrorithm: Polarities for the clickpoints relative to the incoming and outgoing velocity slopes:
        # amps = [1,-1,0,-1,1] # based on five points within the timing gap. 
        # if ((V1 < 0) and (V2 < 0)): # - - starting velocity negative, ending velocity negative
        # amps = [-1,1,0,-1,1]
        # if ((V1 < 0) and (V2 > 0)):    # - + 
        # amps = [-1,1,0,1,-1]
        # if ((V1 > 0) and (V2 > 0)):   # + + 
        # amps = [1,-1,0,1,-1]
        # if ((V1 > 0) and (V2 < 0)):  # + -
        # amps = [1, -1,0,-1, 1]     
        # Check Amean, it seems to be calculating incorrectly. -drb 1/7/2025
        
def generate_points3(T1, A1, V1, T2, A2, V2, Numpoints, Tmed, Tmin, Alow, Ahigh,Amean, cpmean):
    points = []
    # T1 represents beginning time of segment
     # Calculate total time gap representing the length of time that needs to be filled for this segment
    Tgap = T2 - T1
    # Scenarios determine the starting and ending direction of the point amplitudes depending on intial velocities 
#
#
#  Make the amps multipliers velocity-aware. There are however, some glitches at the beginning and ending of the generated miniseed files files. 
#
#    amps = [-1*V1/np.abs(V1),1*V1/np.abs(V1),-1*V2/np.abs(V2),1*V2/np.abs(V2)] # Set velocity direction based on bracketing clickpoints
    amps = [-1*V1/np.abs(V1),1*V1/np.abs(V1),0,1*V2/np.abs(V2),-1*V2/np.abs(V2)] # Set velocity direction based on bracketing clickpoints
#    amps = [0,-1*V1/np.abs(V1),0,-1*V2/np.abs(V2),0] # Set velocity direction based on bracketing clickpoints starting at velocity zero?
    print(f'beginning velocity: {V1} : {amps} : Ending velocity: {V2}')
    
    npoints = 5  # This defines how many points are created.
    clickpoint_time = []
    clickpoint_amplitude = []
    for i in range (0,npoints): # Fix the length to four points for now
        clickpoint_time.append( T1 + (i+1)*(Tgap/(npoints+1))) # Time1 + 1/6th of the timing gap, up to 5/6th of the way to T2
        #
        # The Offset needs to be derived as a function of time, between the mean of the samples BEFORE the gap to the mean of the samples AFTER the gap.
        #
        

        offset =  (A1+((A2-A1))/((npoints+1)/(i+1))) # Break the amplitude into five intervals that begin at A1 + 1/6th difference and end at A2 - 1/6th the difference 
        print(f'A1:{A1} A2:{A2} offset:{offset} for point {i}')
        #clickpoint_amp = amps[i]*(np.random.uniform(Ahigh/2, Ahigh))
#        clickpoint_amp = amps[i]*(np.random.uniform(0, Ahigh))
#        clickpoint_amp = (amps[i]*np.abs(np.random.uniform(Alow, Ahigh)) - Amean) # - offset)
        clickpoint_amp = (amps[i]*np.abs(np.random.uniform(-cpmean, cpmean))) #  - Amean + offset)
#        clickpoint_amp = (amps[i]*(np.random.uniform(Alow, Ahigh)) - Amean) # - offset)
#        clickpoint_amp = (Amean+offset + amps[i]*(np.random.uniform(0, Ahigh)))
#        clickpoint_amp = (Amean+offset + amps[i]*np.abs(Ahigh))
#        print(f'clickpoint_amp: {clickpoint_amp} - Amean: {Amean} + Offset: {offset}  = {clickpoint_amp - Amean + offset}')
#        clickpoint_amplitude.append((clickpoint_amp - Amean + offset)) 
#        clickpoint_amplitude.append(Amean+offset + clickpoint_amp)
        print(f'clickpoint_amp: {clickpoint_amp} + Offset: {offset}  = {clickpoint_amp + offset} and Amean = {Amean}')   
        clickpoint_amplitude.append(clickpoint_amp + offset)
        #clickpoint_amplitude.append(offset + amps[i]*np.abs(np.random.uniform(Ahigh/2, Ahigh)))

    # Adjust A by shifting it so mean is zero
    
    points = []
    for i in range(0,len(clickpoint_time)):
        points.append((clickpoint_time[i],clickpoint_amplitude[i] )) # - np.mean(clickpoint_amplitude))) # This is the mean of the gap.
    print(f'Ahigh: {Ahigh} cm Alow: {Alow}cm Amean: {Amean}cm Mean of clickpoint gap: {np.mean(clickpoint_amplitude)} cm')
    return points

#
# def find_segments (Used by fillgaps()):
# Findsegments finds and reports where there are missing segments of clickpoints representing the time gaps at the top of the minute.
#

def find_segments(clicktime,clickamp,PNEdydx):
    #
    # Isolate and report the missing data segments.
    # Notice that it will sometimes grab "extra" stuff that might or might not be the time mark
    # Parse out only segments that land on a minute mark.
    # Note that PNEdydx (instantaneous velocity) is not accurate for this applicaiton,
    # and a new velocity must be calculated from the points before and after the bracketing clickpoints
    #    segments:
    #    [0] segmentclicktime
    #    [1] segmentinterval 
    #    [2] segment starting clickpoint amplitude
    #    [3] velocity leading into the segment
    #    [4] segment ending clickpoint amplitude
    #    [5] velocity exiting the segment

    Clickinterval = [0]
    segment = []
    for i in range(1,len(clicktime)):
        interval = clicktime[i]-clicktime[i-1] # Clicktime time gaps 
        Clickinterval.append(interval)
        # Only take gaps larger than (slightly less than) 1 second that land on an even 60-second interval 
        # Cut it off at 2.5 seconds as anything longer is likely a signal dropout or an undigitizable section of the waveform.
        # These bracketing clickpoints fo NOT represent inflection points but points where the trace disappears from view.
        # i represents the index for the ending clickpoint of the segment
        # If trapped, report back i-1 as the index for the beginnign of the segment
        # Calculate the matching velocity of the clickpoint. 
        # For velocity through clickpoint[i-1], use (clickamp[i-1] - clickamp[i-2]]) / (clicktime[i-1] - clicktime[i-2])
        # For velocity through clickpoint[i], use (clickamp[i+1]-clickamp[i])/(clicktime[i+1] - clicktime[i]) 
        #
        #                       Valid segment falls on a period of 60 seconds from beginning of waveform, +- 0.25 seconds.
        if (interval) > 0.92 and (interval) < 2.5 and (int(round(clicktime[i-1],0)/60) == round(clicktime[i-1],0)/60) :  # Some intervals are slightly less than one second.
            segment.append([clicktime[i-1],interval,clickamp[i-1],\
           (clickamp[i-1] - clickamp[i-2]) / (clicktime[i-1] - clicktime[i-2]),\
            clickamp[i],(clickamp[i+1]-clickamp[i])/(clicktime[i+1] - clicktime[i]),i-1])        
    return(Clickinterval,segment)

#
# def fillgaps (used by main() and PNE2Clicktime())
#  Find and fill timing mark gaps with a synthesized data segment that looks like neutral data and miminizes artifacts 
# when applying filters and instrument correction to the digitalized signal
#

def fillgaps(clicktime,clickamp,PNEdydx,metadata):
    #
    #
    # If metadata['fillgaps']: Adjust the clickpoint stream to isolate timing gap segments, then correct by generating synthetic clickpoints.    
    #
    #
#
# Analyze the click points and assemble a list of segments where interval exceeds 1 second, and occurs at the top of the minute.
#
# Inputs include:
# clicktime: A list of click points and their time relative to REFtime
# clickamp: The digitized amplitude of the clickpoint in centimeters
# metadata: A list of common variables related to the signal and program processing options.
# PNEdydx : A list of the instantaneuous velocities of each clickpoint in centimeters / second
#
    print(f"The fillgaps def is being executed.")
    Clickinterval,segments = find_segments(clicktime,clickamp,PNEdydx)
    print(f"The number of segments identified for this waveform is {len(segments)}.")

#    segment[0] = segmentclicktime
#    segment[1]segmentinterval 
#    segment[2]clickamp segment starting clickpoint amplitude
#    segment[3] PNEdydx velocity leading into the segment
#    segment[4] clickamp segment ending clickpoint amplitude
#    segment[5] PNEdydx velocity exiting the segment
#    segment [6] index i of this segment 

    segmenttime = []
    segmentvalue = []
    segmentindex = []
    for segment in segments:
        segmenttime.append(segment[0]) # clicktime (seconds elapsed after the first second of REFtime)
        segmentvalue.append(segment[1])# timing gap length in seconds
        segmentindex.append(segment[6]) # segment index
        print(f"At ({segment[6]}:{segment[0]:0.2f} seconds) missing {segment[1]:0.3f} sec data:")
        print(f"    Starting amplitude for this segment: {segment[2]:0.3f} cm  Ending amplitude:{segment[4]:0.3f} cm" )    
        print(f"    Starting velocity for this segment: {segment[3]:0.3f} cm/sec  Ending velocity:{segment[5]:0.3f} cm/sec \n\n" )
        # Calculate the time of projected beginning and ending inflection point, based on an average amplitude for the past
        # ten clickpoints. segment[4] is the start velocity, segment[6] is ending velocity. segment[7] is index for beginning
        # clickpoint and segment[7]+1 is ending clickpoint
        # 
    # This brackets the timing gap by several seconds and returns the max, min, and peak-peak for the sections 
    # immediately before and after the timing gap.
    # On average, there are "n" clickpoints per second in the waveform.
    # Number of clickpoints needed to synthesize the timing gap should be rounded up to the next integer.
    # Try doubling this number.
    # Try using the np.max clickpoints per second.

    for i in range(0,len(segmentindex)):
        print(segmentindex[i],segmenttime[i],segmentvalue[i])
    print("\n\n")
    #
    #  Build the synthesized clickpoints for each segment by generating a grouping of random points with amplitudes within the bounds of adjacent clickpoints.
    #
    frogdna = []
    for i in range(0,len(segmentindex)):
        offset = -5   # Set an arbitrary number of clickpoints before the gap to fifteen clickpoints.
        theend = 2*np.abs(offset)
        if segmentindex[i] < 5:
            offset = 0
        if len(clicktime) - segmentindex[i] < theend: # set the end of the set to what's left
            theend = (len(clicktime) - segmentindex[i])
        clicktimedifference = []
        for j in range(segmentindex[i]+offset,segmentindex[i]+offset+theend): # Go to ten clickpoints past the gap
            clicktimedifference.append(clicktime[j+1] - clicktime[j])
            print(f"{j}: {clicktime[j]} sec. {clickamp[j]} cm.")
        # Print the max and min amplitude values for this (currently 20) clickpoint series
        print(f'For segment at {segmenttime[i]} seconds, with a timing gap of {segmentvalue[i]} in length:')
        print(f'Clicktime min = {np.min(clicktimedifference)} Ave = {np.mean(clicktimedifference)} median = {np.median(clicktimedifference)} std = {np.std(clicktimedifference)} max = {np.max(clicktimedifference)}')
        #print (f'Max amplitude: {np.max(clickamp[(segmentindex[i]+offset):(segmentindex[i]+offset+theend)])} Min amplitude: {np.min(clickamp[(segmentindex[i]+offset):(segmentindex[i]+offset+theend)])}')
        print (f'Peak-peak: {np.max(clickamp[(segmentindex[i]+offset):(segmentindex[i]+offset+theend)]) - np.min(clickamp[(segmentindex[i]+offset):(segmentindex[i]+offset+theend)])}')
        print (f'Mean amplitude: {np.mean(clickamp[(segmentindex[i]+offset):(segmentindex[i]+offset+theend)])}') # Pass this to the generate_clickpoint3
        #
        # Calculate the rate of clickpoints/second for the 20 samples surrounding the timing gap
        # For experimentation, try doubling the number of clickpoints to see how it affects the displacement
        # and velocity waveforms.
        #
        
        rate = 1/np.min(clicktimedifference) 

        numpoints =  math.floor(rate * segmentvalue[i])+1 # number of points for the timing gap
        # calculate maximum number of points necessary to fill the gap
        # Then double it for increasing frequency content
        # Then calculate the time interval represented by each point.
        # To synthesize the signal, let each point land on the top of a half-sin lobe.
        # If one point, time of a period should be twice the signal gap length.
        # If two points, period = signal gap length + a period of 2x signal gap length
        # If three points, 1/3gap + 1xgap + 2xgap 
        #print(theend,(clicktime[segmentindex[i]+offset+theend]-clicktime[segmentindex[i]+offset]-1))
        #print(f'Average clickpoints per second for this segment: {rate} points/second. \n Create {numpoints} clickpoints/second for this segment.')
        #print(f'Number of seconds between each synthesized clickpoint = {segmentvalue[i]/(numpoints+1)}')
        clickpoint_synthetic = []
    #
    #  Build clickpoints for the segments using a random amplitude value 
    #
        # cpmin, cpmax, when there is a slope in the near-DC displacment can result in overly large representitive signals when the signal level is small.
        # Need to apply a linear regression to the sample clickpoints in order to calculate the peak-to-peak for the window.
        # 
        # Send as cpmin, the mean - (peak/peak / 2) and for cpmax, the mean + peak/peak/2 which represent signals immediately above and below the mean
        # Build a list of peak-to-peaks:
    
    #    cpmin = np.min(clickamp[(segmentindex[i]+offset):(segmentindex[i]+offset+theend)]) # Alow : overshoot a little bit
    #    cpmax = np.max(clickamp[(segmentindex[i]+offset):(segmentindex[i]+offset+theend)]) # Ahigh: overshoot a little bit
    #    cpeak2peak = np.max(clickamp[(segmentindex[i]+offset):(segmentindex[i]+offset+theend)]) - np.min(clickamp[(segmentindex[i]+offset):(segmentindex[i]+offset+theend)])
        Amean = np.mean(clickamp[(segmentindex[i]+offset):(segmentindex[i]+offset+theend)]) # the mean of the segment
        A1 = np.mean(clickamp[(segmentindex[i]+offset):segmentindex[i]])
        A2 = np.mean(clickamp[(segmentindex[i]):(segmentindex[i]+offset+theend)])
        cpeak2peak = [] 
        for j in range((segmentindex[i]+offset+1),(segmentindex[i]+offset+theend)):  # Attempt to generate a list of the peak2peak values from clickpoint to clickpoint
            cpeak2peak.append(np.abs(Amean - clickamp[j])) # difference between the clickpoint and the mean of the samples surrounding the gap
            print(j,clickamp[j],np.abs(Amean - clickamp[j]))
        cpmax = np.max(cpeak2peak)
        cpmin = np.min(cpeak2peak)     
        cpmean = np.mean(cpeak2peak)
        print(f'cpeak2peak = {cpeak2peak} , cpmin = {cpmin}, cpmean = {cpmean} cpmax = {cpmax}')
        
        tmin = np.min(clicktimedifference) 
        tmed = np.median(clicktimedifference)
        # rate goes as 7th variable
        #
        #    segment[0] = segmentclicktime
        #    segment[1]segmentinterval 
        #    segment[2]clickamp segment starting clickpoint amplitude
        #    segment[3] PNEdydx velocity leading into the segment
        #    segment[4] clickamp segment ending clickpoint amplitude
        #    segment[5] PNEdydx velocity exiting the segment
        #    segment [6] index i of this segment 
        # generate_points inputs:
        # T1 = time of beginning point      (segment[0]) : segmentclicktime
        # A1 = Amplitude of leading series of clickpoints
        # V1 = Velocity at beginning point  (segment[3]) 
        # T2 = Time of endin point
        # A2 = Amplitude at ending series of clickpoints 
        # V2 = Velocity at ending point
        # rate of points per second
        # tmed: median time gap
        # tmin: minimum time gap
        # cpmin: lowest variation from mean of adjacent clickpoints
        # cpmax: highest variation from mean of clickpoints
        # Amean: The mean of adjacent clickpoints
        # print(segments[i][0],segments[i][1])
        #
        # generate_points3 = Build four or five ( or n) clickpoints, with amplitudes that move from segments[2] to segments[4], but
        #                    with amplitudes of 90% peak from cpmax to cpmin. For instance, with 5 clickpoints:
        #                    Clockpoint[0] occurs at T1+(T2-T1/6) * 1 with amplitude: (segments[3]-segments[2]/6)* 1
        #                    Clickpoint[1] occurs at T1+(T2-T1/6) * 2 with amplitude: (cpmax-cpmin)/2+(segments[3]-segments[2]/6) * 2
        #                    Clickpoint[2] occurs at T1+(T2-T1/6) * 3 with amplitude: -(cpmax-cpmin)/2+(segments[3]-segments[2]/6) * 3
        #                    Clickpoint[3] occurs at T1+(T2-T1/6) * 4 with amplitude: (cpmax-cpmin)/2+(segments[3]-segments[2]/6) * 4
        #                    Clickpoint[4] occurs at T1+(T2-T1/6) * 5 with amplitude: -(segments[3]-segments[2]/6)* 5
        #              # generate_points3(T1, A1, V1, T2, A2, V2, Numpoints, Tmed, Tmin, Alow, Ahigh,Amean)                   
        #points = generate_points3(segments[i][0], segments[i][2], segments[i][3], segments[i][0]+segments[i][1], \
        #                segments[i][4], segments[i][5], rate, tmed, tmin, cpmin, cpmax,np.mean(clickamp[(segmentindex[i]+offset):(segmentindex[i]+offset+theend)]))
        points = generate_points3(segments[i][0], segments[i][2], segments[i][3], segments[i][0]+segments[i][1], \
                        segments[i][4], segments[i][5], rate, tmed, tmin, cpmin, cpmax, Amean, cpmean) # Make the mean halfway between the two clickpoints
        print(f'Points generated: \n{points}')
        frogdna.append(points)   # Add the newly created clickpoints to the list for integration into the main clickpoint file for this component at the timing gaps.
        print('\n\n')
        
    # At this point, we have four parallel lists:
    # segmentindex[] 
    # segmenttime[]
    # segmentvalue[]
    # frogdna[]
    # Insert the new clickpoints into the original clickpoint list.
    # insert them, starting at the end gap and working backward to the first timing gap to preserve the index.
    # With a segmentindex and accompanying synthetic clickpoints, 
    # insert the click points into the clickpoint stream 
    #
    # The segment index contains the clickpoint after which we need to append our new synthetic frog dna
    #

    clicktime_modified = clicktime.copy()  # Make some copies for modification
    clickamp_modified = clickamp.copy()
    print(len(clicktime_modified),len(clickamp_modified))  
    segmentlist = []
    #
    # Reverse the clickpoint fills because we'll be filling the LAST timing gap first, and the first timing gap last.
    #
    for i in range (0,len(segmentindex)):
        segmentlist.append([segmentindex[i],segmenttime[i],segmentvalue[i],frogdna[i]])
    segmentlist.reverse()
    for segment in segmentlist: # Going through each segment, insert it into the original click point list to fill in those gaps 
        Frogdna = segment[3]  # generated frogdna clickpoints from generate_points3
        Frogdna.reverse() # Invert the segment and pack them into the list
        #
        # point by point, end to beginning , insert the list of clickpoint times and amplitudes into the main clicktime/clickamp list 
        # at the index specified at segment[0]+1. This grows the original clicktime/clickamp lists by appending the lists with the synthesized datapoints.
        #
        print('FrogDNA segments being inserted into the clickpoint file:')
        for time,amplitude in Frogdna:
            print(time,amplitude)
            clicktime_modified.insert(segment[0]+1,time) 
            clickamp_modified.insert(segment[0]+1,amplitude)
    #
    # Troubles at the beginning and the end of the file when velocities are high...
    # Set clickpoint[0] and clickpoint[len(ckickpoint)-1] to value of points adjacent to force zero velocity at the beginnign and end of the list.
    #
#    clickamp_modified[0] = clickamp_modified[1]
#    clickamp_modified[len(clickamp_modified)-1] = clickamp_modified[len(clickamp_modified)-2] 
   
    return (clicktime_modified,clickamp_modified)


#
# Optimize function checks for unrealistic spikes and optimizes clickpoint placement to minimize high-velocity
# spikes above the assumed seismometer channel operating bandwidth
# then provides a modified list of clickpoints along with diagnostic lists
# such as the clickpoints adjusted and related amplitudes,
# the testdata derivative, and the spikes that exceed the threshold.
# Only adjust clickpoints that exceed 1/4 the minmax for the 20 adjacent clickpoints as small changes with asymmetry are
# likely real.
# clicktime, clickamp represent click point file for generation of signal.
# streamtime, streamdata = actual obspy stream from the original clickpoint file, 

def Optimize(clicktime,clickamp,streamtime,streamdata,metadata,time_params):
    # This is a hard rule for now, but may change.

        
    tr = Maketrace(streamdata,metadata,time_params) # makes a trace from streamdata, a time-history output from PCHIP
    tr2 = tr.copy()
    tr2.differentiate(edge_order=2) # Create the differential representing Dy/Dt for this should be roughly symmetrical without huge deviations from sample to sample
    tr2.filter("highpass",freq = metadata['breakpoint'], corners = 2) # single-pole filter above operating frequency to find big spikes
    testdata = [(sample) for sample in tr2.data] # convert differential channel to an absolute value as we look for excursions
    spiketime,spikeamp = Find_spikes((testdata),streamtime,float(metadata['threshold'])) # test against a threshold
    #
    # Match requires the spiketime match the clicktime within +- delta
    #
    # find minimum time between click points
    #
    clickminimum = find_minimum_click_interval(clicktime)
    print(f'minimum time between clicks is {clickminimum:0.3f} seconds.')

    Clicktime_to_adjust,Clicktime_to_adjust_index = Match(spiketime,clicktime) # find problem clickpoints that made the spikes
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
        PNEdydx = [] # This one is going to be a crude estimate of direction and slope by going through the list, setting slope
        # as change in amplitude as a change in time between clickpoints. Last clickpoint will use slope of zero.
        metadata = {'comment':'','station':'','network':'','component':'','reftime':'01_JAN_1970_00:00:00.000', \
        'starttime':'01_JAN_1970_00:00:00.000','cf':float(1.0),'tc':float(0.0), 'location':'', 'samplerate':float(100.0),\
        'nsamples':int(0),'polarity':1.0,'threshold':25.0,'breakpoint':16.0,'optimize':'False','shiftlimit':0.085,\
        'fillgaps':'False', 'Nograph':'False'}
        for row in list:
            if len(row) > 0:          # skip blank lines
                r = row[0].split()    # Split up the un-blank lines into discrete values then check them.
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
 #   metadata['samplerate'] = float(metadata['samplerate'])
    metadata['nsamples'] = int(((clicktime[len(clicktime)-1]-clicktime[0])*metadata['samplerate'])+1) # This is the number of samples to generate with interpolation algorithm
    metadata['optimize'] = str2bool(metadata['optimize'])
    metadata['fillgaps'] = str2bool(metadata['fillgaps'])
    first_sampletime    = float(clicktime[0])# Number of seconds ahead of the ref_time in which the first sample actually occurs
    Time_params = Get_time_params(metadata,first_sampletime)
    #
    #                       Delta is the average sample period. It is calculated from the offset time of last sample 
    #                       minus offset of first sample / total number of samples.
    # 
    metadata['delta']     = (float(clicktime[len(clicktime)-1])-first_sampletime)/(metadata['nsamples'])
    print(f"delta has been set to {metadata['delta']} seconds, which is {1/metadata['delta']} Hz")
    metadata['samplerate'] = 1/(float(clicktime[len(clicktime)-1])-first_sampletime)/(metadata['nsamples'])
    errtime = False # There is no error time in a mouse click file.
    #
    # Calculate PNEdydx
    #
    PNEdydx.append(0) # Append a slope of zero for the first data point
    for i in range(1,len(clicktime)): # Go to the second-to-last point 
        dydx = (clickamp[i]-clickamp[i-1]) / (clicktime[i]- clicktime[i-1]) # centimeters / seconds
        PNEdydx.append(dydx)
    PNEdydx.append(0) # Append a slope of zero for the last data point
    return(clicktime,clickamp,PNEdydx,errtime,metadata,Time_params)

#
#  def Export_clickfile: Write the mouse click file representing the vectorized clickpoints of the seismogram to disk.
#
def Export_clickfile(metadata,clicktime,clickamp):
    clickfile = os.path.join(metadata['outfolder'],(metadata['filename']+"_clickpoints.csv"))
    with open(clickfile, mode='w', newline='\n') as mouseclick:
        mouseclick_writer = csv.writer(mouseclick, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        mouseclick_writer.writerow(['COMMENT',metadata['comment']])
        mouseclick_writer.writerow(['STATION',metadata['station']])
        mouseclick_writer.writerow(['NETWORK',metadata['network']])
        mouseclick_writer.writerow(['LOCATION',metadata['location']])
        mouseclick_writer.writerow(['COMPONENT',metadata['component']])
        mouseclick_writer.writerow(['CF',metadata['cf']])
        mouseclick_writer.writerow(['TC',metadata['tc']])
        mouseclick_writer.writerow(['REFTIME',metadata['reftime']])
        mouseclick_writer.writerow(['STARTTIME',metadata['starttime']])
        mouseclick_writer.writerow(['OPTIMIZE',metadata['optimize']])
        mouseclick_writer.writerow(['THRESHOLD',metadata['threshold']])
        mouseclick_writer.writerow(['SHIFTLIMIT',metadata['shiftlimit']])
        mouseclick_writer.writerow(['BREAKPOINT',metadata['breakpoint']])
        mouseclick_writer.writerow(['FILLGAPS',metadata['fillgaps']])
        mouseclick_writer.writerow([])
        mouseclick_writer.writerow(['Clickpoint','Relative_time(seconds)','Amplitude(centimeters)'])
        for i, mouseclicks in enumerate(clicktime,0):
            mouseclick_writer.writerow([i,f'{(clicktime[i]-clicktime[0]):003.003f}',clickamp[i]])
    print ("Mouseclick file successfully written: {0}".format(clickfile))
    return()
#
# Load a wavetrack output file with header.
# The header contains the necessary metadata for creating a miniseed file from the wavetrack project file.
# The stack contains two element list, of time in seconds since first sample, amplitude in centimeters.
#    
def Import_Wavetrack(infile):
    with open(infile,'r') as fin:
        list = csv.reader(fin)
        PNE = []
        metadata = {'comment':'','station':'','network':'','component':'','reftime':'01_JAN_1970_00:00:00.000', \
        'starttime':'01_JAN_1970_00:00:00.000','cf':float(1.0),'tc':float(0.0), 'location':'', 'samplerate':float(100.0),\
        'nsamples':int(0),'polarity':1.0,'threshold':50.0,'optimize':'False','breakpoint':16.0,'filename':'none',\
        'shiftlimit':0.085,'outfolder':'C:\\','interpolation':'pchip','fillgaps':'False','Nograph':'False'}
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
#    metadata['samplerate'] = float(metadata['samplerate'])
    metadata['cf']         = float(metadata['cf'])
    metadata['tc']         = float(metadata['tc'])
    metadata['shiftlimit'] = float(metadata['shiftlimit']) # shiftlimit is listed in seconds.
    metadata['optimize'] = str2bool(metadata['optimize'])   # Boolean value based in interpretation of text field
    metadata['starttime'] = metadata['starttime'].upper()
    metadata['reftime'] = metadata['reftime'].upper()
    metadata['fillgaps'] = str2bool(metadata['fillgaps'])

    #
    #                       Build the time parameters dictionary
    #

    first_sampletime    = float(PNE[0][0])# Number of seconds ahead of the ref_time in which the first sample actually occurs 
    Time_params = Get_time_params(metadata,first_sampletime)

    #                       In New digitizations, first_sampletime is expected to be zero.
    #                       Delta is the average sample period. It is calculated from the offset time of last sample 
    #                       minus offset of first sample / total number of samples.

    Delta     = (float(PNE[len(PNE)-1][0])-first_sampletime)/(len(PNE)-1)
    metadata['delta'] = Delta                # This is the calculated sample interval, which should coincide with the stated sample rate
    metadata['samplerate'] = 1/Delta
    metadata['nsamples'] = len(PNE)          # remember to adjust this value when trimming the file to size
    
    return(metadata,PNE,Time_params)



#
#
# Get Time Parameters: (used by import_wavetrack and import_clickfile) Parse out the start time and reference time, generate UTCDateTime instances
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

#    PNE2Clicktime: Used by main() to extract the click points from the wavetrack output file as well as create
#    the PNE_time/PNE_amp lists representing the digitalized signal resulting from those clickpoints.
#                       Load Data array from the Wavetrack output file
#                       Use PCHIP to provide a list of click times, and interpolated waveform from those click times.

def PNE2Clicktime(metadata,PNE,Time_params):
    # version 20240810 - Also pass PNEdydx outside the def for use in fillgaps() if necessary

    PNEtime = [] # represents time track of the selected dataset starting from the selected start time
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
    #                Process file for click points, then interpolate into smoother waveform with an interpolation method such as pchip

    metadata['nsamples'] = len(PNEtime)    # Re-adjust the metadata to account for the trimming of the excess samples
    PNEamp   = np.arange(len(PNEamplitude),dtype=np.float32)
    for n in range(len(PNEamplitude)):      
        PNEamp[n] = np.float32(PNEamplitude[n]) 
    PNEamp = polynomial(PNEamp, order = 1, plot = False) # remove the offset with a polynomial linear trend .

    # Get the clickpoints
    # You need to bring out the estimated instantaneous derivative at the clickpoint as well.
    # Call it clickdydx
    #
    # Extrapolate the click points and amplitudes from the wavetrack linearly interpolated points.
    # Where the wavetrack changes slope, there's a clickpoint. Determine exactly where the two line segments intercept and calculate the time
    # and amplitude of that line segment intercept point, and make it the clickpoint.
    #
    PNEclicktime,PNEclickamp,PNEdydx,errtime = WT_to_clickpoint(PNEtime,PNEamp)
    #
    # It is HERE that one checks for filling in the timing gaps via metadata parameter fillgaps(clicktime,clickamp,GMTClicktime,PNEdydx,metadata)
    #
   
    if metadata['fillgaps'] == True:
        print(f"fillgaps option is set to {metadata['fillgaps']}:\n")
        PNEclicktime, PNEclickamp = fillgaps(PNEclicktime,PNEclickamp,PNEdydx,metadata)
    #
    #
    
    if metadata['interpolation'] == 'pchip':
        PNE_time,PNE_amp = Pchip(PNEclicktime,PNEclickamp,len(PNEtime))
    if metadata['interpolation'] == 'cubicspline':
        PNE_time,PNE_amp = Cubicspline(PNEclicktime,PNEclickamp,len(PNEtime))
    if metadata['interpolation'] == 'cubichermitespline':
        PNE_time,PNE_amp = CubicHermiteSpline(PNEclicktime,PNEclickamp,PNEdydx,len(PNEtime))
    if metadata['interpolation'] == 'none':
        PNE_time,PNE_amp = nointerpolation(PNEtime,PNEamp)

    #
    #               Prepare a list of UTC corrected times that correspond to the click points.
    # PNEclicktime is a list of clickpoints in UTCDateTime
    PNE_clicktime = [Time_params['St_time']+(sample - Time_params['Starting_sample'])+float(metadata['tc']) for sample in PNEclicktime] 
    Clicktime = []
    for click in PNE_clicktime:
        Clicktime.append(str(click))
    return(PNE_time,PNE_amp,PNEclicktime,PNEclickamp,Clicktime,errtime,PNEdydx)

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

# Line intercept is used to find the exact time of the clickpoint based on where the two slopes intersect.
#  Here is where the following code is sourced and subsequently adapted to PNE2SAC:
#  https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines


def line_intercept(point): # version 20250210
                            # Points are configured as [time(x), amplitude(y), slope(m)]
    error = 0
    point1 = point[0]
    point2 = point[1]
    b1 = (point1[1] - point1[2]*point1[0])  #b = amplitude(y) - slope(m) * time(x)
    b2 = (point2[1] - point2[2]*point2[0])  # this is b: the y intercept.
    m1 = point1[2] #  slope at point # 1 
    m2 = point2[2] #  slope at point # 2
    if m1 == m2:
        error = -1
        interpolated_time = point1[0]+(point2[0]-point1[0])/2.0
        interpolated_amplitude = point1[1]+(point2[1]-point1[1])/2.0
        return (interpolated_time,interpolated_amplitude,0,error)
    x = (b2 - b1) / (m1 - m2) # Calculate the x position of the interection of the lines
    y = m1 * x + b1 # 
    slope = m1 * (x - point1[0])/(point2[0]-point1[0]) + m2 * (point2[0] - x)/(point2[0]-point1[0])
    # in reality, the slope should be weighted depending on the "closeness" between the two points.
    # point1[0] is time at first point, point2[0] is time at second point, x is time of clickpoint between these two.
    # If the calculated point falls more than .010 seconds outside of the wavetrack points, assume
    # there is a discontinuity within wavetrack data and interpolate.
    
    if not point1[0] <= x <= point2[0] and ((x - point1[0] < -.010) or (x - point2[0] > 0.010)) :
        error = x
        x = point1[0]+(point2[0]-point1[0])/2.0 # interpolate a midpoint
        y = point1[1]+(point2[1]-point1[1])/2.0
  #      slope = 0  # uncalculatable slope
    return (x,y,slope,error) # Pass coordinates of the clickpoint plus the miscalculated value for error logging.




def WT_to_clickpoint(time,amp): # input the time and amplitude series. method describes the desired interpolation method
    # version 20250131
    # method = 0 [default] = Pchip interpolation methodology
    # method = 1 = spline but this is not yet implemented.
    # Slope = difference in amplitude across two adjacent samples divided by time difference between the time of those two samples
    slopes = []
    bslopes = []
    fslopes = []
    deltat = [] # [time[1]-time[0]]
    dither_ratio = 6 # Detection threshold for flagging samples that bracket a change in slope representing a clickpoint.
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
            if (np.abs(slopes[i+1] - slopes[i])) > dt/dither_ratio : # dt = samplerate . dt/dither_ratio determines how much 'wiggle' in the slope from sample to sample is allowed before triggering
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
    PNEdydx = [] # Let's attempt to assemble a list of average slopes as a centimeter / second value
    PNEclicktime = []
    PNEclickamp = []
    errtime = []
    for point in points:
        if not point[1]:       # a second point does not exist for the first and last point within the list.
            PNEdydx.append(point[0][2]) # The slope associated with the first or last clickpoint
            PNEclicktime.append(point[0][0])
            PNEclickamp.append(point[0][1])
        else:                  # find the intercept point from the points list and add it to the final list.
            x,y,slope,error = line_intercept(point)
            PNEdydx.append(slope)
            PNEclicktime.append(x)
            PNEclickamp.append(y)
            if error != 0:
                log_error(point,x,y,error) # Process the error
                errtime.append(x) # Record whenever a linerly interpolated point exists.

    return(PNEclicktime,PNEclickamp,PNEdydx,errtime)
#
# def Pchip (used for converting the clickpoint list of vectorized clickpoints into a digital signal)
# (PCHIP stands for Piecewise Cubic Hermite Interpolating Polynomial)
# clues taken from here: https://stackoverflow.com/questions/31221444/intuitive-interpolation-between-unevenly-spaced-points
#
def Pchip(clicktime,clickamp,nsamples): # Bring in clickpoints and output a PCHIPped time series
    # At this point, clicktime represents the time of each click point and clickamp represents amplitude of each mouse click
    x_data = np.array(clicktime)
    y_data = np.array(clickamp)
    Pchip_time = np.linspace(min(x_data), max(x_data), nsamples) # From the oldest time to the youngest time, space out 'n' points, as specified in len(slopes)
    bi = interpolate.PchipInterpolator(x_data, y_data)
    Pchip_amplitude = bi(Pchip_time)
    return(Pchip_time,Pchip_amplitude)

def Cubicspline(clicktime,clickamp,nsamples): # Bring in clickpoints and output a PCHIPped time series
    # At this point, clicktime represents the time of each click point and clickamp represents amplitude of each mouse click
    x_data = np.array(clicktime)
    y_data = np.array(clickamp)
    Pchip_time = np.linspace(min(x_data), max(x_data), nsamples) # From the oldest time to the youngest time, space out 'n' points, as specified in len(slopes)
    bi = interpolate.PchipInterpolator(x_data, y_data)
    bi = interpolate.CubicSpline(x_data, y_data)
    Pchip_amplitude = bi(Pchip_time) # This is the time axis onto which the interpolated points are mapped. Returns the amplitudes on this axis
    return(Pchip_time,Pchip_amplitude)

def CubicHermiteSpline(clicktime,clickamp,clickdydx,nsamples):
    # At this point, clicktime represents the time of each click point and clickamp represents amplitude of each mouse click
    x_data = np.array(clicktime)
    y_data = np.array(clickamp)
    dydx = np.array(clickdydx)
    Pchip_time = np.linspace(min(x_data), max(x_data), nsamples) # From the oldest time to the youngest time, space out 'n' points, as specified in len(slopes)
    # Before we can use this, we must generate an array of matching derivatives (dydx) for each click point
    bi = interpolate.CubicHermiteSpline(x_data, y_data, dydx)
    Pchip_amplitude = bi(Pchip_time) # This is the time axis onto which the interpolated points are mapped. Returns the amplitudes on this axis
    return(Pchip_time,Pchip_amplitude)

def nointerpolation(time,amp):
    # Create a stream of linearly interpolated points straight from the original wavetrack output file.
    return(time,amp)

#
# Convert the digitalized list of clickpoints into a SAC stream for export.
#
def Create_sacstream(b,metadata,time_params,outfolder):
    #                                          Create the SAC stream
    St_time = time_params['St_time']  + float(metadata['tc'])          # Shortcut for filename below
#    Delta     = metadata['delta'] # (1/float(metadata['samplerate']))
    Network   = metadata['network']
    Component = metadata['component']
    Location  = metadata['location']
    Stname    = metadata['station']
    t         = SACTrace(data = b)         
                                 # set the SAC header values
    t.scale  = 1.0               # Set the scale for use with DIMAS software
    t.delta  = metadata['delta']
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
    print (" File successfully written: {0} with a start time of {1} and sample period of {2:0.3f} seconds".format(outfile,St_time,t.delta))    
    sacfile.close()
    return(t) #  This is the trace that was just written to disk.

#
# Convert the SAC trace into a miniseed trace for eventual export. Why do it this way? Well, this WAS once a PNE2SAC code, after all.
#
#
def Make_obspy_trace(tr,metadata): # This converts the sac trace tr into an obspy trace then returns an obspy trace
    tr2 = tr.to_obspy_trace()
    #tr.stats.comment    = metadata['comment'] 
    tr2.stats.network   = metadata['network'] 
    tr2.stats.station   = metadata['station'][:7] 
    tr2.stats.channel   = metadata['component'][:3]
    tr2.stats.location  = metadata['location']
    tr2.stats.starttime = starttime
    tr2.stats.delta = metadata['delta'] # 1/float(metadata['samplerate'])
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
    tr.stats['delta']= metadata['delta']
    tr.stats['network'] = metadata['network']
    tr.stats['station'] = metadata['station'][:7]
    tr.stats['location'] = metadata['location']
    tr.stats['channel'] = metadata['component'][:3]
    tr.stats['starttime'] = time_params['St_time']+float(metadata['tc']) # Add the time correction here
    return(tr)

#
# plot_optimized makes a plot with all the information, from non-optimized/optimized comparisons
# to clickpoints, to some stats about the digitization.
# It also calculates some statistics about the streams.
#
def plot_optimized(metadata,time_params,outfolder,PNE,tr1,tr2,errtime,streamtime,clicktime,Clicktime_to_adjust,spiketime,clickamp,spikeamp,testdata,newtestdata,C2amp):
    PNEtime = []
    PNEamplitude = []
    offset = False

    for data in PNE:
        if float(data[0]) >= (time_params['Starting_sample']):  # Starting_sample:  # Discard the unwanted samples from beginning of file 
            if not offset:
            #    print(f'Setting offset to {offset}')
                offset = float(data[1])
             #   print(f'Setting offset to {offset} while data[1] = {data[1]}')
            PNEtime.append(float(data[0])) # This is the relative time
            PNEamplitude.append(-1*(float(data[1])-offset)) # This is amplitude in centimeters

    St_time = time_params['St_time']
    tr4 = tr2.copy()
    tr4.differentiate(edge_order=2)
    tr4.filter("highpass",freq = metadata['breakpoint'], corners = 2) # single-pole filter above operating frequency to find big spikes
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
    if not errtime:
        outputselect = (f"Final output file represents the {waveselect} waveform.")
    else:
        outputselect = (f"Due to existing discontinuitie(s), no output files are generated.")

    fig, axs = plt.subplots(3)
    fig.suptitle(f"{optimized}\n{shiftlimit}\n{errors}\n{adjusteds}\n{outputselect}")
    fig.text(0.06, 0.5, 'Amplitude (centimeters)', ha='center', va='center', rotation='vertical')
    fig.set_size_inches(15.0,10.0)
    xmin = 0
    xmax = clicktime[len(clicktime)-1]
    
#                                       For zoomed-in plots that emphasize the clickpoints in need of adjustment, set the boundaries thus:

    if len(Clicktime_to_adjust) > 1:
            if Clicktime_to_adjust[0] > 30:
                xmin = Clicktime_to_adjust[0]-30
            if  (clicktime[len(clicktime)-1] - Clicktime_to_adjust[len(Clicktime_to_adjust)-1])  > 60: 
                xmax = Clicktime_to_adjust[len(Clicktime_to_adjust)-1]+60
    else:
        xmin = clicktime[0]
        xmax = clicktime[len(clicktime)-1]
        
    #
    # Set the limits for each of the three subplots
    #
    axs[0].set_xlim([xmin,xmax]) # Set the upper plot to zoom in ten seconds before the first adjusted clickpoint to ten seconds after the last adjusted clickpoint
    axs[1].set_xlim([clicktime[0],clicktime[len(clicktime)-1]]) # Set the middle plot to show the whole waveform
    axs[2].set_xlim([xmin,xmax]) # Set the threshold plot to the same as the top plot.
    #
    # Annotate the three subplots
    #
    axs[0].annotate('Wavetrack points (black), Original PCHIP(blue), optimized waveform(red). Discontinuities in green.',xy=(xmin+(xmax-xmin)/4.0, np.min(tr2.data)))
    axs[1].annotate('optimized waveform (red)',xy=(xmin+(xmax-xmin)/2.0, np.min(tr2.data)))
    axs[2].annotate('Trigger channel before(blue) vs. after(red)',xy=(xmin+(xmax-xmin)/2.0, np.min(testdata)))
    #
    # Add an error bar onto the top plot for any discontinuities that might be present.
    #
    axs[0].bar(errtime, height = np.max(tr1.data)*2, width=.2, bottom=np.min(tr1.data)*2, align='center',color = "green",fill=True)
    #
    # Plot the original PNE stream of interpolated points, if they exist.
    #
    if(len(PNEtime)):
        axs[0].scatter(PNEtime,PNEamplitude,color='black',s=2)
    #
    # Plot the data
    #
    axs[0].plot(streamtime,tr1.data,label = "Original trace",color='blue',linewidth = 2) # original trace in blue
    axs[0].scatter(clicktime,clickamp,color = 'blue',s=8) # original clickpoints
    axs[0].plot(streamtime,tr2.data,color='red',label = "Modified trace",linewidth = 1) # optimized trace overlay
    axs[1].plot(streamtime,tr2.data,color='red') # optimized trace n second plot
    if C2amp: # C2amp represents the clickpoints that is are candidates to adjust, if the optimizer is turned on.
        axs[0].scatter(Clicktime_to_adjust,C2amp,color = 'red',s=24) # adjusted clickpoints
    if spiketime:
        axs[2].scatter(spiketime,spikeamp,color = 'red',s=8) # Where threshold is exceeded
    if testdata: # trigger channel for finding the clickpoints that exceed the threshold set by the operator.
        axs[2].plot(streamtime,testdata,color='blue', linewidth = 2) 
    axs[2].plot(streamtime,newtestdata,color='red',linewidth = 1)
    #
    # Add documentation to the plot
    #
    ftime = (f"{St_time.year}.{St_time.month:02.0f}.{St_time.day:02.0f}.{St_time.hour:02.0f}.{St_time.minute:02.0f}.{St_time.second:02.0f}")
    filename = metadata['network']+"."+metadata['station']+"."+metadata['location']+"."+metadata['component']+"."+ftime
    filenametitle = (f"Digitized channel name: {os.path.join(outfolder,filename)}.mseed\n{metadata['comment']}\n PNE2SAC version {__version__}")
    fig.text(0.5, .0625, filenametitle, ha='center', va='center' , rotation='horizontal')
    plt.savefig(os.path.join(outfolder,filename+".png"))
    if metadata['Nograph']: # For batch processing of export files the plot is saved to disk and not pushed to the user's screen.
        plt.close()
    else:
        plt.show()
    return()
    
    
def main():
#                               Parse the command line switches and set the defaults
    Folder = False  # Treat the input as a single file unless the -folder command line switch is present.
    SAC = True      # Always make sac files in this version.
    MSEED = True    # Always generate a miniseed file, unless -display is present as a command line switch
    QCcheck = True  # Always perform a QCcheck and do not output files if QC fails.
    Sps = 100       # PCHIP interpolated output will be 100 sps
    outputfile_defined = False
    filelist = []
    Wavetrack = True
    Clickfile = False
    Fillgaps = False
    Nograph = False
    Logfile = True   # Changed for Dan's batch processing , normally is False
    infile = ''
    infolder="./"                              # Assume user CWD is the same folder as the target file until specified otherwise
    outfolder = ".\\"                         
    extension = '.txt'
    
    args=parser.parse_args()                  # parse the arguments supplied to the program. 
    if args.clickfile:                        # By default the extension is .txt so change it for clickfiles to .csv   
            Clickfile = True                  # Switch to clickfile processing and don't try to extract clickpoints.
            Wavetrack = False
            extension = ".csv"                # change the default extenstion from .txt to .csv 
    if extension in args.filename:            # Options are clickfile, wavetrack output file, or a folder. 
            filelist.append(args.filename)    # a single file has been speficied within the filelist field.
            infolder = args.filename[:args.filename.rfind("\\")]  # The input folder for this file.
            infile = args.filename            # The file name of the input file.
    if args.folder:                           # If -folder option, then the target is not a file name but is instead infolder.
            try:
                Folder = True                 # This means the first instance following folder represents the folder name.
                infolder = args.filename
                filelist = os.listdir(args.filename)
            except:
                print("Warning! No folder specified within arguments. Check syntax.")
    if args.nograph:                          # Don't output the graph at the end (for batch processing of wavetrack files)
            Nograph = True
    if args.logfile:                          # Set text output to a log file within the destination folder.
            Logfile = True
    if args.optimize:                         # Enable optimization algorithm to fine tune clickpoints
            metadata['optimize'] = True       # Superceded if optimize is specified within the file header
    if args.fillgaps:
            metadata['fillgaps'] = True       # this will be superceded if fillgaps is specified within the file header.
    if args.display:                          # Only display the waveforms, don't overwrite the existing sac and miniseed files within the folder.
            SAC = False
            MSEED = False

    Sps = float(args.sps)                     # Set the pchip sample rate which by default is 100 sps (about 10X oversampling) for good amplitude resolution in the time domain.

    if args.filename.rfind("\\") > 1:                 # unless a folder is in the name.
        outfolder = args.filename[:args.filename.rfind("\\")]

    for n in range(len(filelist)):            # filelist normally has but one entry except in the case of being a folder in which case the code will concurrently convert files.
        if extension in filelist[n]:          # Ensure that the file is a valid input file with valid extension, (.txt or .csv depending on the switch settings)
            if len(filelist)>1:
                infile = os.path.join(infolder,filelist[n])
            outfolder = ".\\"                          # Assume user has specified a local file
            if infile.rfind("\\") > 1:                 # unless a folder is embedded within the name.
                outfolder = infile[:infile.rfind("\\")]
#
# Establish link to an standard output file for the log.
#
            if (Logfile or Nograph):                                # Turn off the graph at the end so that the code can be run in batch mode from Jupyter Notebook
                print(f" Outfolder set to: {outfolder} \nRedirecting all output to {infile[:infile.rfind('.')]}.log")
                sys.stdout = open(infile[:infile.rfind(".")]+'.log','w')
#
#            Document the PNE2SAC version and run time for this code execution within the log file.
#
            print(f'\nPNE2SAC Version {__version__} executed on {time.asctime(time.gmtime())} UTC\n')
#
#                 Import the wavetrack output file or the clickfile
#
            if Wavetrack:
                metadata,PNE,time_params = Import_Wavetrack(infile)
                if Fillgaps:                           # Set the fillgaps switch True
                    metadata['fillgaps'] = True
#
#               Process wavetrack points into a time-history stream at the specified sample rate
#                    
                streamtime,streamamp,clicktime,clickamp,GMTClicktime,errtime,PNEdydx     = PNE2Clicktime(metadata,PNE,time_params) # Extract the clickpoints from the WT output file
                streamdata        = np.arange(len(streamamp),dtype=np.float32)
                for n in range(len(streamamp)):        #   Load the array with time-history data
                    streamdata[n] = np.float32(streamamp[n])
                St_time = time_params['St_time'] + float(metadata['tc'])
                print(f"Importing Wavetrack file {infile}:")
#
#               Import the clickfile of previously generated clickpoints. 
#
            if Clickfile: 
                PNE = []
                clicktime,clickamp,PNEdydx,errtime,metadata,time_params = Import_clickfile(infile)
                if Fillgaps:
                    metadata['fillgaps'] = True
                #
                # Here, we check for if timing gaps should be filled via fillgaps(clicktime,clickamp,GMTClicktime,PNEdydx,metadata)
                # Return a modified version of clicktime,clickamp, metadata, if fillgaps = True
                #
                if metadata['fillgaps']:
                    clicktime, clickamp = fillgaps(clicktime,clickamp,PNEdydx,metadata)
#
#               Process the clickpoints into a time-history stream at the specified sample rate.
#               
                streamtime,streamamp     = Pchip(clicktime,clickamp,metadata['nsamples'])
                streamdata        = np.arange(len(streamamp),dtype=np.float32) # b = time history intended for placing within the miniseed stream
                for n in range(len(streamamp)):                                #   Load the array with time-history data
                    streamdata[n] = np.float32(streamamp[n])
                St_time = time_params['St_time']+ float(metadata['tc'])
                print(f"Importing Clickfile {infile}:")
            #
            #                             Graphs are turned on/off depending on whether this is run as a stand-alone code of as part of a list of many files.
            #
            if Nograph:
                metadata['Nograph'] = True
            else:
                metadata['Nograph'] = False
#
#           Print out diagnostics for the conversion of this file
#
            print( "\nSample rate = {0:2.3f} \nNumber of samples = {1}".format( \
                  1/float(metadata['delta']),len(streamdata)))
            if len(PNE):
                print( "First_sampletime = {0} seconds, last_sampletime = {1:2.3f} seconds.".format( \
                       float(PNE[0][0]),float(PNE[len(PNE)-1][0])))
                print(f"Starting sample for the time slice occurs at {float(time_params['St_time']-time_params['Rf_time'])} seconds.") 
            print(f"Time correction = {metadata['tc']} seconds.")
            print(f"starttime = {time_params['starttime']}")
            #
            # errors within the clickpoint calculations due to corrupted wavetrack clickpoints (typically two clickpoints on the same pixel on time axis)
            #
            if errtime:            # Discontinuities within the file mean there are errors in the wavetrack project file in need of correction. Dont output faulty files.
                QCcheck = False
                print(f"There are {len(errtime)} discontinuities in this file.\n")
                if len(errtime):
                    print("Repair discontinuities in Wavetrack output file, then re-run PNE2SAC.")
            #
            #               Calculate optimized trace even if it is unused.
            #               Calculate the optimized version of the clickpoint stream 
            #               
            new_clicktime,clickamp,Clicktime_to_adjust,C2amp,testdata,spiketime,spikeamp = \
            Optimize(clicktime,clickamp,streamtime,streamdata,metadata,time_params)
            New_streamtime,New_streamamp     = Pchip(new_clicktime,clickamp,metadata['nsamples']) # Make new PCHIP interpolation
            tr1 = Maketrace(streamdata,metadata,time_params) # makes a trace from old streamdata, a time-history output from PCHIP
            tr2 = Maketrace(New_streamamp,metadata,time_params) # generate a new trace from revised clickpoints.
            tr3 = tr1.copy()
            tr3.filter("highpass",freq = metadata['breakpoint'], corners = 1) # single-pole filter above operating frequency to find big spikes
            tr4 = tr2.copy()
            tr4.differentiate(edge_order=2)
            tr4.filter("highpass",freq = metadata['breakpoint'], corners = 2) # single-pole filter above operating frequency to find big spikes
            newtestdata = [sample for sample in tr4.data] # Newly corrected stream
            #
            # At this point, we have the original trace tr1, optimized trace tr2, and the trigger channel from the original at tr3.
            #
            # Generate the file name for this dataset and place it in metadata to pass it to other defs for use in saving and plotting.
            #
            ftime = (f"{St_time.year}.{St_time.month:02.0f}.{St_time.day:02.0f}.{St_time.hour:02.0f}.{St_time.minute:02.0f}.{St_time.second:02.0f}")
            filename = metadata['network']+"."+metadata['station']+"."+metadata['location']+"."+metadata['component']+"."+ftime
            metadata['filename'] = filename
            metadata['outfolder'] = outfolder
            print(f"\nNumber of triggers in original channel = {len(spiketime)}; Number of clickpoints adjusted = {len(Clicktime_to_adjust)}")
            print(f"Maximum allowable shift of a clickpoint is defined as {metadata['shiftlimit']} seconds.")
            print(f"mean of original stream trigger = {np.mean(np.abs(testdata)):.2e} mean of modified stream trigger = {np.mean(np.abs(newtestdata)):.2e}")
            print(f"max spike amplitude of original stream = {np.max(np.abs(testdata)):.2e} max spike amplitude (modified) stream = {np.max(np.abs(newtestdata)):.2e}")
            print(f"%reduction in spike amplitude: {(1-np.max(np.abs(newtestdata))/np.max(np.abs(testdata)))*100.:02.1f} %\n")
            outfile = ''
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
                        print(f"Using the optimized trace to create miniseed file. \nSample period = {tr.stats['delta']}")       
                    else:
                        tr = tr1 
                        print(f"Using the original trace to create miniseed file. \nSample period = {tr.stats['delta']}")
                    tr.write(outfile,format = "MSEED")
                    print(f"File succesfully written: {outfile}\n\n")
            else: # Flagged QCcheck has failed so no output is generated.
                print('Quality check failure; Check for discontinuities and re-edit the Wavetrack file.')
                print('No SAC or miniseed files written to disk for this waveform.')
            #
            # Plot the resulting streams, including the original wavetrack output, error bars, clickpoint list, and interpolated waveform, plus spikes, discontinuities, etc.
            #
            plot_optimized(metadata,time_params,outfolder,PNE,tr1,tr2,errtime,streamtime,clicktime,Clicktime_to_adjust,spiketime,clickamp,spikeamp,testdata,newtestdata,C2amp)
            #
            # Restore the screen output so script can terminate
            #
            if (Logfile or Nograph):
                sys.stdout.close()
                sys.stdout = sys.__stdout__ 

#
# Check and run the main function here:
#
if __name__ == '__main__':
  main()