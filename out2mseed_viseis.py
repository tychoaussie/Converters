# -*- coding: utf-8 -*-

'''The MIT License (MIT)

Copyright (c) 2013 Daniel Burk
during my time on campus during a summer of delightful albeit, uncompensated labor.
Michigan State University.

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
__version__ = "20180725"
__license__ = "MIT"

# -*- coding: utf-8 -*-
# 20180725 - Version customized for use with VISEIS data as found in early station Talya
# 20180523 - Lots more troubleshooting and dealing with timing problems. Re-written the way the timing is calculated so more
#            accurate measurements are being taken. Studied how the timing gets out of whack- Finally gave up in trying to 
#            correct for unknown timing and now I quarantine any file that has grossly miscalculated in terms of where the
#            timing mark is, since I cannot tell if it's advanced by a second or retarded. Evidence has been found where timing
#            marks will sometimes spontaneously stop for most of a second, or double-strike a few hundred milliseconds from
#            one another. This is discriminated by analyzing the 'range' of samples found between the timing marks. Samples 
#            between timing marks should always be within one sample of another if sample rate is consistent and timing is
#            consistent. Any files where the range exceeds 2 samples means the timing mark has jumped forward or backward, or
#            samples have been added/dropped. These files are quarantined.
#            I think GPS might have been jumping online and off-line, sending timing responsibility to the computer which simulated
#            timing pulses or something. So if the file doesn't have a correctly calculated timing mark, I now discard it by moving the
#            file into a quarantine subfolder. Additionally, the 'dayfile' option is now available, which copies the data into a
#            dayfile subfolder where big files are broken into 24-hour chunks. It remains to be seen if it is useful however,
#            since the discarding of files with suspect timing now seems to make a 24 hour continuous file pretty much improbable.
# 20180510 - Add channel awareness and make up a unique temp folder name to prevent old data from rolling into the merge
# 20180406 - Lots of troubleshooting to account for time instabilities, inconsistent channel naming conventions,
#            etcetera. CWB didnt like the new merged data files so they are explicitly being designated as STEIM2 MSEED files.
#            Program is being retrofitted to output into log files that accompany each converted folder.
#            It also started getting utilities to ensure the maximum file size is one 24 hour day file.
#
#            Note that in time correction, it is highly unusual to have to time shift more than twice in a row.
#            Re-engineer the code and if a file needs a time shift, keep track of it's name in a list.
#            Think of some way to keep a "perma-shift" in the data from happening.
#
#            Bug: Sometimes the program will shift time backward one second and all subsequent files will shift.
#                This time shift occurs after a file has been 'discarded', and previous valid stop time exists.
#                The start time is assumed to be correct. Start times are always suspect when the timing pulse falls
#                within the first 200 milliseconds. So, discard the file if frevious start time is zero AND frac_second
#                is less than 200 milliseconds. 
#                New discriminator has been added: Data[0][2] represents the timing pulse seconds remainder from the digitizer.
#                If the file is bad, chances are this will not synchronize to the calculated timing pulse Frac_second.
#                Thus a new boolean metric called timing (True/False) is created based on timingpulse being within 50 msec of the
#                calculated fractional second.
#
#
# 20180330 - Integrate the miniseed merge function that requires use of local drive for miniseed
#            buffer. This speeds up the conversion to about 30 minutes per CD. Timng looks okay now
#            sample command to run it on a single file:
#            out2mseed ./SEY03N42 RY ../../sdcms/SEY/SEY03N42 > ./SEY03N42/SEY03N42_conversion.log
#            The local drive is designated as c:/seismo/temp and you must manually clear it in-between runs.
#
#
# 20180329 - Bring in file logging, scan and avoid zero file length OUT files. Catch and isolate
#            OUT files that have timing problems.
# 20170519 - This is a linux version, with several improvements to the timing validation.
#            Ver. 0519 also corrects for any mildly corrupt vdaq.txt files.
#
# 20170517 - Removed the need for SAC file writing from the code and now just go straight to miniseed.
#
# 20170516 - Stop trying to calculate sample rate and just trust the one procided in the vdaq
#             file. THe sample rate GPS pulse is not always consistent and sometimes has glitches
#             causing the cample rate calculator to fail. There's no easy fix, so just go with
#             the published sample rate, which seems to be more or less correct.

# 20170510 - Corrected a timing bug in the calculation of sample rate and seconds remainder.
# OUT files incorrectly report start time when the timing pulse occurs within the first four
# or five samples. This requires knowledge of the previous file's calculateed end time.
# As a result, this program must have a time continuous list of files. 
 

# 20170203 - update to use the new SAC file engine from LANL and OBSPY verson 1.01
# Now outputs to Miniseed with STEIM2 compression, and creates a more sensible, trackable
# file name
# 20150817 - updated to run on latest Obspy version plus bring in channel names into
# the SAC filename. Also correct the start time bug by using the encoded start time
# from the first two samples of the timing signal.

import sys, os, glob, csv, time, datetime, string, subprocess, numpy as np
from obspy.core import read, Trace, Stream, UTCDateTime

# from obspy.io.sac import SACTrace


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


#class OUTconvertError(Exception):
#
#   def __init__(self, error_name, error_text):
#       self.args = (error_name, error_text)
#       self.error_name = error_name
#       self.error_text = error_text
#
#   def __str__(self):
#       return repr('{0}: {1}'.format(self.error_name, self.error_text)


class OUTconvert(object):
    '''OUT2MSEED is a utility for converting OUT binary files as generated by
       the Russian seismic networks that use Symmetric Research seismic hardware.
       The OUT files are converted to miniseed format for integration within other
       seismic analysis tools. At present, only OUT files generated in version 1.14 
       and later are compatable with the converter. Future iterations of this code
       should also include early OUT files as well as conventional Symres DAT files.

       Syntax:     OUT2SAC Targetfolder Networkcode (Destinationfolder)[crawl] [dayfile]
       
       where:      Targetfolder is the folder containing the OUT file(s).
                    Networkcode is the two letter designator for your seismic network.
                    Destinationfolder is where you want the miniseed files to be written.
       [Optional]   crawl will go through Targetfolder looking for subfolders of OUT files.
       [Optional]   dayfile will take the miniseed files and split them into 24-hour long day files.
       
       It is a requirement that the target folder contains the file "vdaq.txt" 
       that describes the data file header information.

       Useage requirements: You must have installed both Python27, ObsPy, 
       and NumPy on your machine in order to run this package. Dependencies
       include the os, csv, time, string, numpy, and obspy modules.

       Typical useage:
       c:\Python27\scripts> python OUT2mseed.py c:/data/mystation/2018/datafolder LM ./mystation/seed
       c:\Python27\scripts> python OUT2mseed.py c:/data/mystation/2018 LM ./mystation/seed crawl dayfile

       '''
#
#    OUT2mseed is a utility for converting ascii files generated by Symmetric Research 
#       and Russian network .OUT files into a ObsPy compatable format. 
#       Once reformatted, ObsPy may be used to convert the data into SAC or even
#       Miniseed format.
#       
#       Useage: Activate Python's ObsPy package, then from the command prompt
#       specify the target file and the output filename. If no output filename
#       is provided, the default will be called output.ASC in your current
#       working directory. If a directory is specified, all OUT files found within
#       will be converted.  The target directory in all cases must contain at least
#       one OUT file and one vdaq.txt file with appropriate header information.#
#
#       Syntax: OUT2MSEED input_directory NETDesignator output_directory [crawl] #
#
#       crawl is an optional command line switch that assumes input directory contains
#       subfolders that actually contain the data.
#
#       Typical useage:
#       <ObsPy> C:\Python27\scripts> python OUT2mseed.py /../infil /../outfolder crawl
#       <Obspy> /home/public/SeismicData/sdc/sdcraw/SEY/> python out2mseed.py ./ RY ../../../sdcms/SEY/ crawl
#      python ../../../out2mseed.py ./ RY w:/SeismicData/sdc/sdcms/SEY/SEY03N39/ > ./SEY03N39_conversion.log


#
# Dayfile: input the file, figure out the dates, and then write all of the day files into the dayfile subfolder.
#          Dayfile needs to be modified to pass through multiple time segments if they reside in the same day.
#          The object of this code is to split the huge multi-week files that merge sometimes makes
#          into 24 hour-long files. It is not fully implemented yet.
#
def dayfile(folder,infile):
    success = False
    try:
        dayfolder = os.path.join(folder,"dayfile")
        if not os.path.exists(dayfolder):
            os.mkdir(dayfolder)
        st = read(os.path.join(folder,infile))
        start= st[0].stats.starttime.date
        end = st[0].stats.endtime.date 
        date = start
        for i in range((end-start).days+1): 
            filetime = str(date)+"T"+str(datetime.time(0))
            dt = UTCDateTime(date)
            a = st.slice(dt,dt+(3600*24))
            filename = "{0:04d}_{1:02d}_{2:02d}_{3:02d}_{4:02d}_{5:02d}_{6}_{7}_{8}_{9}.mseed".format( \
                a[0].stats.starttime.year,a[0].stats.starttime.month,a[0].stats.starttime.day,a[0].stats.starttime.hour, \
                a[0].stats.starttime.minute,a[0].stats.starttime.second,a[0].stats.network,a[0].stats.station, \
                a[0].stats.location,a[0].stats.channel)
            a.write(os.path.join(dayfolder,filename),format="MSEED",encoding="STEIM2")
            date += datetime.timedelta(days=1)
            success = True 
            print "Writing day file: {}".format(os.path.join(dayfolder,filename))
    except:
        print "Dayfile operation unsuccesful for {}".format(os.path.join(folder,infile))
    return(success)

#
# Merge: Merge the thousands of three-minute files that were placed in the temporary folder 
# and place the consolidated miniseed files into the destination folder. These files are
# separated by any data discontinuities (gaps) that might exist.
#
def filemerge(infolder,outfolder,Channel):
    success = False
    print "Number of channels to merge = {}".format(len(Channel))
    if len(Channel) > 0:
        print "Merging the files from {0} into {1}\n".format(infolder,outfolder)
        if not os.path.exists(outfolder):
            os.mkdir(outfolder)
        extension = ".mseed"
        channels = Channel # ["BHN","BHE","BHZ","TIM*"] # Consider passing channel names into the Definiton
        print "channels imported for merging: {}".format(channels)
        for channel in channels:
            if not (channel == 'UNK'):    
                st = read(os.path.join(infolder,"*"+channel+".mseed")) # This operation takes some time.
                print "\nThere are {0} streams in channel {1}\n".format(len(st),channel)
                #   Open the 'gaps' file that the converter creates when it finds a gap more than hundred milliseconds
                #   Compliant mseed files will read YYYY_MM_DD.
                gaps = []                    # Find and index the meaningful gaps
                gap = st.get_gaps()
                for i in range(0,len(gap)):
                    if np.abs(gap[i][6]) > 0.100: # Arbitrarily set the allowable gap max length to 100 milliseconds. 
                        print "Meaningful gaps at:{0} for time of {1}".format(i,gap[i])
                        gaps.append(gap[i])

                beginning_index = 0         # Initialize the stream index
                ending_index = 0      # If there are no gaps, it's all one big roll
                print "\n"
                for gap in gaps:
                               # 5th element is the gap ending time that coincides with the stream start time following the gap
                    for i in range(beginning_index,len(st)):
                        if (st[i].stats.starttime == gap[5]):  # Check the stream start time against the gap ending time
                            ending_index = i                   # Ending stream is exactly one less than this number
                            print "st[{0}]stats.starttime is a match to gap: {1}".format(ending_index,gap[5])
                    print "Index = st[{0}] = end of continuous stream with start time of {1}".format(ending_index-1,st[ending_index-1].stats.starttime)
                    print "Merging streams from {0} through {1}".format(beginning_index,ending_index-1)
                    a = st[beginning_index:ending_index].merge(fill_value='latest')
                    filename = "{0:04d}_{1:02d}_{2:02d}_{3:02d}_{4:02d}_{5:02d}_{6}_{7}_{8}_{9}.mseed".format( \
                                a[0].stats.starttime.year,a[0].stats.starttime.month,a[0].stats.starttime.day,a[0].stats.starttime.hour, \
                                a[0].stats.starttime.minute,a[0].stats.starttime.second,a[0].stats.network,a[0].stats.station, \
                                a[0].stats.location,a[0].stats.channel)
                    a.write(os.path.join(outfolder,filename), format = "MSEED", encoding = "STEIM2")#, reclen = 512)
                    print "Successfully written file {0}\n".format(os.path.join(outfolder,filename))
                    beginning_index = ending_index

                print "Merging last stream from {0} through {1}".format(ending_index,len(st))

                a = st[ending_index:len(st)].merge(fill_value='latest')
                filename = "{0:04d}_{1:02d}_{2:02d}_{3:02d}_{4:02d}_{5:02d}_{6}_{7}_{8}_{9}.mseed".format( \
                            a[0].stats.starttime.year,a[0].stats.starttime.month,a[0].stats.starttime.day,a[0].stats.starttime.hour, \
                            a[0].stats.starttime.minute,a[0].stats.starttime.second,a[0].stats.network,a[0].stats.station, \
                            a[0].stats.location,a[0].stats.channel)
                a.write(os.path.join(outfolder,filename), format = "MSEED", encoding = "STEIM2")#, reclen = 512)
                print "Successfully written file {0}\n".format(os.path.join(outfolder,filename))
                success = True
    return(success)
#
#
# Returns a list of lists of out files, separated when outfiles are bad with zero bytes.
#
#
def outlist(indir,extension): 
    folderlist = sorted(os.listdir(indir)) #Get the folder and sort it by name
    outfiles = []
    subfiles = []
    print "In Outlist, the number of incoming files = {0}".format(len(folderlist))    
    for file in range(0,len(folderlist)): # Retrieve the filesize information for each folder
        
        if extension in str.lower(folderlist[file]): # If the file is an .OUT file, build an index of filenames and file sizes
            filesize = os.path.getsize(os.path.join(indir,folderlist[file]))
            if filesize < 512:  
                print "File {0} is being quarantined for being a short file with a size of {1}".format(folderlist[file],filesize)
                quarantine(folderlist[file],indir) # Put it in quarantine as a bad file
                if len(subfiles):
                    outfiles.append(subfiles)
                subfiles = [] #  After dumping the list to outfiles, re-initialize it.
            else:
                subfiles.append(folderlist[file])
    if len(subfiles):  # Finish off the list
        outfiles.append(subfiles)
            
    print "There are {0} lists of files within the outfiles list".format(len(outfiles))
    for i in range(0,len(outfiles)):
        print "List item {0} contains {1} files".format(i,len(outfiles[i]))
    return(outfiles) # Returns a list of lists of non-zero length OUT files, assumed to be continuous.


def quarantine(infile,indir):  # File has been deemed bad, with weird dates within and thus needs to be manually examined.
    qfolder = os.path.join(indir,"quarantine")
    if not os.path.exists(qfolder):
        os.mkdir(qfolder)
    os.rename(os.path.join(indir,os.path.basename(infile)),os.path.join(qfolder,os.path.basename(infile)))
    print "Moving {0} to {1}".format(infile,os.path.join(qfolder,os.path.basename(infile)))
    return()

def mkgmtime(t):         # This function makes an inverse for the time tuple to convert back to a floating point epoch
    """Convert UTC tuple to seconds since Epoch"""
    return time.mktime(t)-time.timezone

def SubDirPath(infolder):  # Enable the program to find a subdirectory of folders that presumably contain OUT files.
    return filter(os.path.isdir, [os.path.join(infolder,f) for f in os.listdir(infolder)])  

def vdaq(headerfile):
    #
    #vdaq seeks the headerfile and loads the appropriate variables as defined by the Magadan
    # network vdaq.txt, a modification of symmetric research's dat file format.
    
    with open(headerfile,'r') as head:
        headdata = csv.reader(head)
        header = {}
        for row in headdata:
            r = row[0].split()
            # Load a dictionary of header file information
            try:
                header[r[0]] = str.upper(r[1])
            except:
                print "Warning: {0} has some unfilled data fields.".format(headerfile)
    return header

def Process_Folder(outfolder,filelist,extension,indir,S_time,S_time_Frac_second,Network):
    Channel = [] # Default to zero length list in case process_folder doesn't find any files
    if not os.path.exists(outfolder):
        print "Creating the folder {}".format(outfolder)
        os.mkdir(outfolder)

    filecount = 0
    Outlist = outlist(indir,extension) # If no files are found, the UNK list of channels will be returned by this def.
    for i in range(0,len(Outlist)): # Run through the parsed filelists that avoid the zero length files
        Tdifference = []                # Note that outlist file names already contain the full path
        for n in range(0,len(Outlist[i])):
            if extension in str.lower(Outlist[i][n]): # filelist[n]):
                Prev_E_time = S_time
                Prev_E_time_Frac_second = S_time_Frac_second

                if len(Outlist[i]) == 1: # Assume the list is discontinuous if there's only one file in there.

                    Prev_E_time = time.gmtime(0.000) # Assume the previous start time is invalid due to discontinuity.
                    Prev_E_time_Frac_second = 0.0
#               You need to add the incoming directory 'indir' to the file name when crawling
                infile = os.path.join(indir,Outlist[i][n])
#                infile = Outlist[i][n]
                print "\nConverting file {0} of [{1}:{2} of {3}]: {4}".format(filecount,i,n,len(Outlist[i]),infile)
#                   Convert returns the start time for the next anticipated continuous file 
                S_time,S_time_Frac_second,Tdifference,Channel = convert(Prev_E_time,Prev_E_time_Frac_second,Tdifference,infile,outfolder,Network)  # Convert the out file into sac files
                filecount = filecount+1
    return(Channel)

#
#      CONVERT:
#      Soubroutine that opens the binary file infile, uses the header information generated 
#      by the vdaq subroutine, and creates an ObsPy compatable ASC file.
# 

def convert(Prev_E_time,Prev_E_time_Frac_second,Tdifference,infile,outfolder,Network):

            #
            #       Look for the header file by parsing out the directory. 
            #       If it doesn't exist, the program will error.
            #       Pass in the previous file's end time to compare to the start time of this file
            #       Bring in the infile, and the specified out file, and the network name.
            #

#    print "The infile specified into the conversion is {}".format(infile)
#    print "The outfile specified into the conversion is {}".format(outfile)
    folder = os.path.split(infile)[0]
#    outfolder = infile[:infile.rfind('/')+1]
    headerfile = os.path.join(folder,'VISEIS.TXT')
    hexfile = infile   #  
    header = vdaq(headerfile)
    Samplecount = header['RecordPts:']
    Channel = []
    for i in range(0,8):
        Channel.append(header["Ch{}ID:".format(i)])


        # Extract the start time to the nearest second from the file name
        # File name is an established standard of 14 characters
        # hexfile[-18:-4] represents st.tiome to nearest second
        # 20130314000509
        # Note!! This is computer system time, NOT the start time. So don't do this.
        # Use the start time as encoded within the timing channel (channel 0) as found
        # within the first two samples.

        #    St_time = time.strptime(hexfile[-18:-4],"%Y%m%d%H%M%S")
 
            # Import the binary data
            # Each channel sample comprises of four bytes
            # Epoch time is taken from bytes 1, 0, 13, 12 in that order.
            # Create a data type object of eight channels each of which consist of a 32bit integer

    dt = np.dtype([(header['Ch0Number:'],np.int32),(header['Ch1Number:'],np.int32),(header['Ch2Number:'],np.int32),(header['Ch3Number:'],np.int32) \
                  ,(header['Ch4Number:'],np.int32),(header['Ch5Number:'],np.int32),(header['Ch6Number:'],np.int32),(header['Ch7Number:'],np.int32)])

            # Calculate sample rate and seconds remainder for adding onto file start time.
            # Load timing signal into an array and calculate mean. 
            # Find first sample representing the first positive edge trigger that's greater than sample 5
            # Note that if signal starts high, it must drop low before counting.
            # Count the number of excursions where timing signal goes initially high, starting with the second timing signal
            # and en
            # Find the first sample where gps tick-mark goes high.
            # If tickmark is already high on the 4th sample, go to the next tick mark and count back.


    Data = np.fromfile(hexfile,dtype = dt)      # load all data from the binary file using our specified format dt
    try:
            # Data[0][0] represents MSBaLSBa 0000 of epoch start time from gps
            # Data[1][0] represents MSBbLSBb 0000 of epoch start time from gps
            # Epoch start time must be arranged thus: MSBa LSBa MSBb LSBb            
            # Note that rest of the Data is simply a list of raw counts from the ADC and is a 32 bit integer value.
            # It would be nice if each channel was converted to a measurement in terms of volts as a float value.
            # The Symres PAR4CH system is 24 bits for +-10V, and the USB4CH is listed as +-4V
            # but practical measurements in the lab put it more like +-8V. 
            # Therefore, this code is going to ASSUME a nominal value of 0.94 microvolts / count.
            # This converter will convert raw counts into millivolts using this gain factor. Future versions
            # will enable us to input the gain factor as it becomes available.
            #
            #
            # Channelgain = [0.94*1e-6,0.94*1e-6,0.94*1e-6,0.94*1e-6] # volts per count
            #
        Fail = False               # If there's a problem, fail the file conversion and stash it into quarantine. 
        GPS = []                   # declare our GPS stream which will be loaded from the Data
        Latch = False
        Count = -1                                  # First time Count is incremented is on tic mark number zero.
        Initial_sample = 0
        Final_sample = 0
        Frac_second = 0.0
        Sps = 0.0
        units = ['Volts   ','Volts   ','Volts   ','Volts   ','Volts   ','Volts   ','Volts   ','Volts   ']
        comment = ['TIME    ','Velocity','Velocity','Velocity','Velocity','Velocity','Velocity','Velocity']
        data = []
        data.append(Data[0][0]) # MSB of start time from file
        data.append(Data[1][0]) # LSB of start time from file
        data.append(Data[2][0]) # Milliseconds elapsed since the start of the previous second

        for n in range(0,len(Data)):                  # Load the GPS array
            GPS.append((Data[n][0]))
        for n in range (0,4):
            GPS[n] = GPS[4]   # Zero out the first four samples that were overwritten by the encoded timing.       
        GPS = GPS - np.median(GPS)                    # Mean removal
                  # zero out the epoch time tics.   
        Gpsmean15 = (np.nanmax(GPS)/50.0)            # Use a trigger value that is significantly smaller than max impulse

#        print "Length of Data = {0} versus length of GPS = {1}".format(len(Data),len(GPS))
#        for n in range(0,20):
#            print "GPS{0} = {1}".format(n,GPS[n])	
#        print "Latch = {}".format(Latch)
#        print "GPS maximum = {}".format(np.nanmax(GPS))
#        print "GPS median = {}".format(np.median(GPS))
#        print "GPS Mean = {}".format(np.mean(GPS))
#        print "Gpsmean15 = {}".format(Gpsmean15)
        Last_count = 0
        Samples_per_pulse = []              # The number of samples between timing marks should not deviate but by one sample
                                            # Any greater deviation means timing pulses are untrusted.
                                            # Check to see if the signal started out on a high pulse
        if GPS[4] > Gpsmean15:
            Latch = True                    # Set latch as rising edge has been missed

        for n in range (5,(len(GPS))):
            if (Latch == True):             # Look for falling edge to reset latch
                if GPS[n] < Gpsmean15:
                    Latch = False
            else:             
                if GPS[n] > Gpsmean15:      # When the timing signal exceeds the threshold, it is considered the rising edge.
                    Latch = True            # Rising edge found so set latch
                    Count += 1
                    if Last_count > 0:              # and increment edge count starting at zero.
                        Samples_per_pulse.append(n - Last_count) # Keep track of how many samples elapse between falling edges.
                    Last_count = n      # Samples between falling edges should be consistent, else timing is suspect.
                                            
                    if Initial_sample == 0:
                        Initial_sample = n  # Set the first known good rising edge
                    else:    
                        Final_sample = n    # Keep updating last known good rising edge


        CalculatedSps = float(int(float(Final_sample-Initial_sample)/Count)*1000)/1000
        Sps = float(header['SampleRate:'])

#        print 'Sample rate set to {}'.format(Sps) 

# Use the vdaq sample rate as above method turned out to be unreliable when
# the GPS signal was glitchy.

        Delta = float(int((1/Sps)*1000))/1000   # Round the Delta to the nearest millisecond.


            #                       Calculate time remainder which equals 
            #                 1000 milliseconds - (#samples before first impulse)

        if (Initial_sample - Sps) == 4 : # Meaning, there are more samples to the edge than
                                    # can be accounted for in one second

            Frac_second = 1 - ((Initial_sample - Sps-1)/Sps)
            print "Initial Sample - SPS = {0}. Fraction of second calculated to {1} sec.".format((Initial_sample - Sps),Frac_second)

        elif  (Initial_sample - Sps) == 3 or (Initial_sample - Sps) == 2: # Meaning, there are more samples to the edge than
                                                                    # can be accounted for in one second

            Frac_second = 1 - ((Initial_sample - Sps-1)/Sps)
            print "Initial Sample - SPS = {0}. Fraction of second calculated to {1} sec.".format((Initial_sample - Sps),Frac_second)
    #           print "Backing up the timing by one second to compensate for OUT file timing glitch."
    #           data[1] = data[1]-1 # Decrement the timer by one second, as the original digitizer
                            # missed the mark.

        elif (Initial_sample - Sps) == 1: # meaning the onset of the edge occurs on the first sample
            Frac_second = 0.0



        else: 
            Frac_second = 1 - ((Initial_sample-1) / Sps)
    #            print "Initial Sample - SPS = {0}. Frac = {1} sec.".format((Initial_sample - Sps),Frac_second)

    #
    #       Validate the timing impulse. data[2] represents the frac_second as sourced from the file. It should match the calculated timing
    #       Do this validation whenever the end time is set to zero.
    #
    # Here is where we really need to fix timing. timepulse tells us that if the Frac_Second is greater than the apparent
    # measurement from data[2], then an overflow has happened and it is now mistakenly advanced to the next second marker.
    # First, subtract one from the seconds counder, data[1]. Then, add it to the milliseconds counter data[2]/1000.
    # Then, subtract Frac_second. The result is the amount you need to remove from the timing in order to get things right.
    # The result is that we subtract a second and use the calculated Frac_Second as the remainder.

#        if Frac_second < 0.0:
#        Frac_second = 1.0 - Frac_second # For those times when it has rolled over 
 
        timing = True
        file_timestamp = long(int(data[0])<<16|int(data[1]))
    
#        if (data[2]/1000.0)-Frac_second < -0.100:     # If timepulse is negative, need to borrow a second from data[1] 
#            filefrac_second = 1+(data[2]/1000.0)
#            timepulse = filefrac_second - Frac_second
#            data[1] = data[1]-1 # Borrow a second from the data[1]
#            print "(1) fraction of second from file:   {0:1.3f} seconds ".format(filefrac_second)
#            print "(1) fraction of second (calculated):{0:1.3f} seconds".format(Frac_second) 
#        else:
        filefrac_second = data[2]/1000.0
        timepulse = (np.abs(1-filefrac_second) - np.abs(1-Frac_second)) #  File fraction of second - calculated fraction of second
        print "fraction of second from file:   {0:1.3f} seconds ".format(filefrac_second)
        print "fraction of second (calculated):{0:1.3f} seconds".format(Frac_second) 
#        endif right here
        print "Reported file timing is different than calculated timing by {0:1.3f} seconds.".format(timepulse)

        if (np.amax(Samples_per_pulse)-np.amin(Samples_per_pulse)) > 2:
            print "{0} Median Samples between timing pulses: {1} range = {2} samples.".format(len(Samples_per_pulse),np.median(Samples_per_pulse),np.amax(Samples_per_pulse)-np.amin(Samples_per_pulse))
            print " Timing is suspect."
            timing = False # Dont trust this file and quarantine it.
    #
    #
    #

        timestamp      = long(int(data[0])<<16|int(data[1])) # Assemble them into the timestamp
        St_time        = time.gmtime(timestamp) # Convert them into a tuple representing start time to nearest second   
        file_starttime = mkgmtime(time.gmtime(file_timestamp))+(data[2]/1000.)# Initial timing from file
        starttime      = mkgmtime(St_time)+Frac_second # This is the proposed start time of the file in epoch with the fraction
        timepulse      = file_starttime-starttime
        print "Difference between file start time(fst) and calculated start time(cst): (fst - cst) = {0:1.3f}".format(timepulse)
        prev_endtime = mkgmtime(Prev_E_time) + Prev_E_time_Frac_second

    #    print  "New end time = {0}{1}".format(time.strftime("%Y_%m_%d_%H_%M_%S_",time.gmtime(prev_endtime)),Prev_E_time_Frac_second)

    #
    #           If there exists a previous end time, it means we trusted the previous file to have a good start and end time.
    #           Vet this new file against this previous time.
    #

        if prev_endtime > 0.999: #   # If the start time is less than one, it indicates this is the first file of the sequence.
    #
    #                                    Begin processing the files that follow a reference file with an assumed good strt time. 
    #  
            print "Proposed start time = {0}{1:3.0f}".format(time.strftime("%Y_%m_%d_%H_%M_%S_",St_time),Frac_second*1000.)
            time_difference = starttime - prev_endtime
    #            print "Fraction of second from file                        = {0:1.3f} seconds".format(filefrac_second)
    #            print "Calculated Fraction of a second for true start time = {0:1.3f} seconds".format(Frac_second)
    #            print "Difference: a possible indicator of timing quality  = {0:1.3f} seconds".format(timepulse)
            print "Previous end time   = {0}{1:3.0f}.       Difference = {2:1.3f} seconds".format( \
                 time.strftime("%Y_%m_%d_%H_%M_%S_",time.gmtime(prev_endtime)),\
                 Prev_E_time_Frac_second*1000, \
                 time_difference)

            if (0.94 > np.abs(timepulse) > 0.100): # Let the timing marks drift by 100 milliseconds, then quarantine the file
                print "Timepulse ({0:0.3f} sec.) drift exceeds 100 msec tolerance. Timing integrity has been lost. ".format(timepulse)
                Fail = True
            elif (1.970 > time_difference > 0.95): # This means the proposed file is too far forward and a second has to be removed
                St_time = time.gmtime(mkgmtime(St_time)-1)
                print "The start time for this file has been shifted back one second."
            elif (-1.970 < time_difference < -0.95): # The proposed file start time is too far backward and one second must be added
                St_time = time.gmtime(mkgmtime(St_time)+1)
                print "The start time for this file has been shifted forward one second."
            elif (-0.20 < time_difference < 0.0):
                if Frac_second < 0.980:
                    Frac_second = Frac_second+0.020 # add one sample onto the start time
                else:
                    Frac_second = 0.0
                    St_time = St_time + 1

            elif (np.abs(time_difference) >  1.97):    # The proposed file start time is beyond normal limits.
                Fail = True
    #
    #                           Otherwise it is a file that has followed a timing gap and must be evaluated to see if it's got good timing. 
    #                           It's start time is suspect if the first sample falls within the first 200 milliseconds so don't use it.
    #                           Fractions of seconds greater than 200 milliseconds MUST match our calculated fractional second if it's
    #                           to be considered as having 'good' time.
    #
        elif ( 5 > (Initial_sample-Sps) > 0 or np.abs(timepulse) > 0.07): # Timing starts less than 100 milliseconds from the beginning of the second so dont trust it
            print "Initial_sample = {0}, Sps = {1}, Initial_sample - Sps = {2}".format(Initial_sample,Sps,Initial_sample-Sps)
            print "Either Initial Timing mark is within {0:0.3f} seconds of file start without a ref time or".format(Frac_second)
            print "Timing mark difference between encoded timing and measured timing = {0:0.3f} seconds exceeds threshold and timing is untrusted".format(timepulse)
            Fail = True

        if timing == False: # timing has been checked above
            print "The timing is untrusted in this file. Timing mark difference = {0:0.3f} seconds".format(timepulse)
            print " Range of samples between timing marks is {0} samples".format(np.amax(Samples_per_pulse)-np.amin(Samples_per_pulse)) 
            Fail = True
    #
    #                     Calculate the rest of the times
    #

        E_time = time.gmtime(mkgmtime(St_time)+len(Data)/Sps+Frac_second)
        E_time_Frac_second = Frac_second+(len(Data)/Sps-int(len(Data)/Sps))

        if E_time_Frac_second >= 1.000: # Dont let it roll over
            E_time_Frac_second = E_time_Frac_second - 1.0

        Start_time = time.strftime("%Y_%m_%d_%H_%M_%S_",St_time)
        End_time   = time.strftime("%Y_%m_%d_%H_%M_%S_",E_time)
        Start_time +=str.format("{0:0.3f}",Frac_second)[2:]
        End_time +=str.format("{0:0.3f}",E_time_Frac_second)[2:]

        if int(len(Data)) <> int(Samplecount): # THe OUT File terminated prematurely.

            print " Warning! Mismatch. Number of samples found: {0} Number of samples reported: {1}".format(len(Data),Samplecount)
            print "File is being quarantined."
            print " Predicted time discontinuity, so end time is being set to zero."
            E_time_Frac_second = 0.0
            E_time = time.gmtime(0.000)
            Fail = True

        filesec = int(infile[-6:-4])
        filemin = int(infile[-8:-6])
        filehour = int(infile[-10:-8])
        filetime = filehour*3600+filemin*60+filesec

        Timediff = filetime-St_time.tm_sec-St_time.tm_min*60-St_time.tm_hour*3600-Frac_second
        Tdifference.append(Timediff)

    #        if not (1.1 > Timediff/np.mean(Tdifference) > 0.9):
    #            print "Warning! Timing has varied more than 10% from the computer clock timing."
    #            print "Timediff = {0}  Mean = {1} seconds".format(Timediff,np.mean(Tdifferrence))
    #            print "File is being quarantined and will not be written to miniseed."
    #            E_time_Frac_second = 0.0
    #            E_time = time.gmtime(0.000)
    #            Fail = True

        et = len(Data)/Sps # Number of samples divided by the sample rate = number of second worth of time-history
    #        print "Sample rate : {0:0.3f}".format(Sps)
    #        print "Ratio Sps/CalculatedSps = {0:0.3f}".format(Sps/CalculatedSps)

    #        print "Calculated samples per second, based on tic marks and samples = {0:0.3f} s/sec".format(CalculatedSps)
        if not (1.1 > (Sps/CalculatedSps) > 0.9):
            print "Warning! Calculated sample rate is more than 10% different than base rate.( {} )".format(Sps/CalculatedSps)
            print "File is being quarantined and will not be written to miniseed."
            E_time_Frac_second = 0.0
            E_time = time.gmtime(0.000)
            Fail = True

    #        print "Ratio Sps/CalculatedSps = {0:0.3f}".format(Sps/CalculatedSps)
    #        print "Delta: {0:8.6e}".format(Delta)
    #        print "{} seconds worth of time-history in this file.".format(et)
    #        print "Number of samples found: {0} Number of samples reported: {1}".format(len(Data),Samplecount)
        print "Official start_time: {}".format(Start_time)
        print "Calculated End_time: {}".format(End_time)
    #        print "Time difference between encoded file name and calculated start time from file: {} Seconds".format(Timediff)
    #        print "Time stamp = {0} {1}".format(data[0],data[1])
    #        print "First impulse found at sample: {}  VALUE: {}\n\n".format(Initial_sample,GPS[Initial_sample])

    #        print "GPS5 = {}".format(GPS[5])
    #        print "GPS Median 2.5 = {0}".format(Gpsmean15)
    #        print "Final impulse found at sample: {} VALUE: {}".format(Final_sample,GPS[Final_sample])
    #        print "Total samples between tic marks = {}".format((Final_sample-Initial_sample))
    #        print "Fraction of second to add: {} second".format(Frac_second)
    #        print "Total count of tickmarks: {}".format(Count)
        if Fail:
            print "The failure flag for this file has been set. It will not be written to file."

            # At this point, we have our header information in an index, and we have calculated the true sample rate, 
            # We have extracted the true start time and we've 
            # verified the true second remainder for placing into the start time.


    #
    #               Begin writing the channels out as miniseed files
    #
    #
        if not Fail:
        #                      Create each stream
            for i in range(0,8):
                stats = {'network': Network, 'station': header['A-DInfo:'], 'location': '',\
                         'channel': Channel[i], 'npts': len(Data), 'sampling_rate': Sps,\
                         'mseed': {'dataquality': 'D'}} 
                stats['starttime'] = UTCDateTime(mkgmtime(St_time)+Frac_second)

                b = np.arange(len(Data),dtype=np.float32)   #   Establishes the size of the datastream
                for n in range(len(Data)):        #   Load the array with time-history data
                    b[n] = np.float32(Data[n][i]) #   Convert the measurement from counts to volts.

                f2 = os.path.join(outfolder,Start_time+"_"+Network+"_"+header['A-DInfo:']+"_{}.mseed".format(Channel[i]))
   
                if St_time.tm_year <> 1970:
          
                    if Channel[i] !="UNK":   # We do not write streams in which the channel name is UNK      
                        st=Stream([Trace(data=b, header = stats)])
                        st[0].data = st[0].data.astype('int32')  # dtype it for use with steim2
                        st[0].write(f2,format="MSEED", encoding = 'STEIM2')
                        print " File successfully written: {0} with {1} samples.".format(f2,len(st[0].data))
                else:
                    print "There are problems with the date in this file so it is not being written to disk."
                    print "Reported year = {}".format(St_time.tm_year)
                    E_time_Frac_second = 0.0
                    E_time = time.gmtime(0.000)
                    quarantine(infile,folder)
            print"\n\n"
            return (E_time,E_time_Frac_second,Tdifference,Channel) # Pass the calculated end time so that it can be used for validating the next file 
        else:
            print "Timing is suspect within this file so it is not being written to disk.\n End time is being set to zero."
            quarantine(infile,folder)
            Fail = False                       # Reset the flag
            E_time_Frac_second = 0.0           #
            E_time = time.gmtime(0.000)        # Begin the next file with the assumption that it's start time is correct.
            return (E_time,E_time_Frac_second,Tdifference,Channel)

    except:
        print " Warning! File {0} has fatal errors.\n Number of samples within file =  {1}\n".format(hexfile,len(Data))
        print " Predicted time discontinuity, so end time is being set to zero."
        print " File is being relocated to quarantine in {0}".format(folder)
        try:
            quarantine(hexfile,folder)
        except:
            print "Can no longer find file {0} for quarantine operation.".format(hexfile)
        E_time_Frac_second = 0.0
        E_time = time.gmtime(0.000)
        return (E_time,E_time_Frac_second,Tdifference,Channel)

def main():

    # Go process the command line switches
    # If the switch for processing a whole directory is called, first catalog the diretory
    # then interatively run through the conversion program for as long as the list exists.
    # If the switch isn't valid, use the previous CLSs to process one specific file.
    # This code will not iterate through bottom directories. 
    # Actually, perhaps a better method would be to NOT embed the crawler into this program
    # but instead use csh and bash and pipe the outputs into thte converter

    #           MAIN PROGRAM BODY
    #  Parse the command line switches
    optioncount = len(sys.argv)
    outputfolder_defined = False
    Dayfile = False
    filelist = []
    indir=""
    outfolder=os.getcwd()+"/out2mseed_output/"
    extension = '.out'
    S_time = time.gmtime(long(0)) # Initialize it to time of the creation of the universe
    S_time_Frac_second = 0.0
    Network = ""
    crawl = False
    tempfolder = os.path.join("./",'temp_'+time.strftime("%Y%m%d%H%M%S"))

    if optioncount > 1:
            # Parse out the command line arguments
            # First argument = incoming file folder
            # second argument = incoming file folder's two digit NETWORK CODE
            # third argument = Outgoing file folder (optional)
            # Fourth argument = "crawl" - This keyword means the first argument is an incoming file folder
            # that does not consist of OUT files, but instead consists of folders containing OUT files.
            # If invoked, Crawl will build, within the outgoing file folder location, an identical subfolder
            # full of miniseed files.

            # C:> out2mseed c:\seismo\RY_MOMR\2013\2013_02\ RY \\MSU-IRIS1\public\SeismicData\NERSP\RY_YBGR\2013\2013_02\ >>out2mseed_output.log
        if optioncount == 6:
            if 'crawl' in sys.argv[5].lower() or 'crawl' in sys.argv[4].lower():
                crawl = True
            if 'dayfile' in sys.argv[5].lower() or 'dayfile' in sys.argv[4].lower():
                Dayfile = True
            outputfolder_defined = True
            outfolder = sys.argv[3] # this is the location to place all subsequent directories
            Network = sys.argv[2]
            indir = sys.argv[1]#+"/"
            filelist = sorted(os.listdir(indir)) # Collect a list of all files in the directory
            infile = filelist[0]

        if (optioncount == 5):
            if sys.argv[4].lower() == "crawl":
                crawl = True
            if sys.argv[4].lower() == "dayfile":
                Dayfile = True
            outputfolder_defined = True
            outfolder = sys.argv[3] # this is the location to place all subsequent directories
            Network = sys.argv[2]
            indir = sys.argv[1]#+"/"
            filelist = sorted(os.listdir(indir)) # Collect a list of all files in the directory
            infile = filelist[0]

        if optioncount == 4:
            outputfolder_defined = True
            outfolder = sys.argv[3]
            Network = sys.argv[2]
            indir = sys.argv[1]#+"/"
            filelist = sorted(os.listdir(indir)) # Collect a list of all files in the directory
            infile = filelist[0]
                 
        if optioncount == 3:    # out2sac infile/directoryname outfile/directoryname
            if ("/" in sys.argv[2]) or ("/" in sys.argv[2]):          
                outputfolder_defined = True
                outfolder = sys.argv[2]
            else:
                Network = sys.argv[2]
                outfolder = sys.argv[1]#+"/"

            indir=sys.argv[1]#+"/"
            filelist = sorted(os.listdir(indir)) # Collect a list of all files in the directory                
            infile = filelist[0]    # initialize for the first iteration

            print "Infolder and outfolder are specified as '{0}' and '{1}'".format(indir,outfolder)
            
        elif optioncount == 2:     # Out2sac inputfile/targetdirectory 
                                   # use default name for target file and directory
            indir = sys.argv[1]#+"/"
            filelist = sorted(os.listdir(indir)) # Collect a list of all files in the directory
            infile = filelist[0]
            outfolder = indir


                                  # Parse through the directory and convert all out files
                                  # If indir, outdir have not yet been specified,
                                  # Specify them.
        if Network == "": 
            Network = raw_input('Please enter the two letter network designator. ')

        Network = (Network[:2]).upper()  # Assume network was typed in correctly and make it a two letter designator
#        print "Total number of arguments: {}".format(len(sys.argv))
#        print "Dayfile = {}".format(Dayfile)
#        print "crawl = {}".format(crawl)
#        for i in range(0,len(sys.argv)): print i,sys.argv[i]
        if crawl:
            # indir = The folder containing a list of subfolders
            # outfolder = The target folder that will contain the created subfolders of miniseed files
            # insubfolder[] = the folder of out files to process
            # outsubfolder = the folder of created miniseed files
            #
            insubfolder = sorted(SubDirPath(indir))   # Make a list of the folders within the source folder 
            if not os.path.exists(tempfolder): # Create the target folder that will contain the subfolders
                print "Creating tempfolder {}".format(tempfolder)
                os.mkdir(tempfolder)

            for j in range(len(insubfolder)): # process the source subfolders
                filelist = sorted(os.listdir(insubfolder[j]))# Make a list of the files within the source subfolder

                # Strip insubfolder name and glue it to the outsubfolder as a target for miniseed
                temp_subfolder = os.path.join(tempfolder,os.path.dirname(insubfolder[j])[insubfolder[j].rfind('/')+1:])
                if not os.path.exists(temp_subfolder):
                    print "Creating temp_subfolder {}".format(temp_subfolder)
                    os.mkdir(temp_subfolder)
                tempsubfolder = os.path.join(temp_subfolder,os.path.basename(insubfolder[j]))
                if not os.path.exists(tempsubfolder): # Create the target folder that will contain the subfolders
                    print "Creating temp_subfolder {}".format(tempsubfolder)
                    os.mkdir(tempsubfolder)
                    #
                # At this point, we have the source subfolder name, the target subfolder name, 
                # the list of files in the source subfolder. Now process the source subfolder.

                print "There are {0} files listed in {1}".format(len(filelist),insubfolder[j])
                if len(filelist)> 0 : # some subfolders might not contain OUT files.
                    Channel = Process_Folder(tempsubfolder,filelist,extension,insubfolder[j],S_time,S_time_Frac_second,Network)
                    success = filemerge(tempsubfolder,outfolder,Channel)
                    if success:                            # Clean up the temp folder to prepare for next crawl
                        #print "tempsubfolder is designated as: '{}'".format(tempsubfolder)
                        for fl in glob.glob(os.path.join(tempsubfolder,"*.mseed")):
                        #    print "Attempting to remove file {}".format(fl)
                            os.remove(fl)
                        os.rmdir(tempsubfolder)
                        S_time = time.gmtime(long(0)) # Re-Initialize to time of the creation of the universe
                        S_time_Frac_second = 0.0
#                    if Dayfile:# A bug in dayfile: The previous folders keep getting re-converted, so dont use this
                                # option in crawl until it is fixed. You probably want to remove the non-dayfiles
                                # after they have been split.
#                        folderlist = sorted(os.listdir(outfolder))
#                        for file in folderlist:
#                            if '.mseed' in file:
#                                success = dayfile(outfolder,file)        

        else:                            # Standard mode of operation without directory crawling
            print "There are {0} files listed in {1}".format(len(filelist),indir)
            if len(filelist) > 0: # Don't try processing files if there's no files to process
                Channel = Process_Folder(tempfolder,filelist,extension,indir,S_time,S_time_Frac_second,Network)
                success = filemerge(tempfolder,outfolder,Channel)
                                              # Clean up the temp folder
                if success:
                    #print "tempfolder is designated as: '{}'".format(tempfolder)
                    for fl in glob.glob(os.path.join(tempfolder,"*.mseed")):
                    #    print "Attempting to remove file {}".format(fl) 
                        os.remove(fl)
                    os.rmdir(tempfolder)
                if Dayfile: 
                    folderlist = sorted(os.listdir(outfolder))
                    for file in folderlist:
                        if '.mseed' in file:
                            success = dayfile(outfolder,file)
    else:
        print "Useage: OUT2MSEED infolder Network outfolder (optional)crawl"
        print "Or, OUT2MSEED infolder"
        print "No input folder has been specified."
        print "Be sure input folder contains .OUT files and a vdaq.txt station descriptor file."
        print len(sys.argv)

#Call the main loop
#
if __name__ == '__main__':
  main()

# <codecell>


