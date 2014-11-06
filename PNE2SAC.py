# -*- coding: utf-8 -*-

'''The MIT License (MIT)

Copyright (c) 2013 Daniel Burk

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
__version__ = "20141106"
__license__ = "MIT"

# -*- coding: utf-8 -*-

# Today, on November 5th, 2014 we lost a family member to illness. We will
# always remember him.

import os, csv, sys, numpy as np, time, string
from obspy.core import read, Trace, Stream, UTCDateTime
from obspy.sac import SacIO

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


def load(infile):
    with open(infile,'r') as fin:
        list = csv.reader(fin)
        rowcnt = 0
        stack = []
        header = []
        for row in list:
            r = row[0].split()
            if rowcnt < 7:
                header.append(r)
                rowcnt +=1
            else:
                stack.append(r)
                rowcnt+=1
                
        return(header,stack)




def main():
#           MAIN PROGRAM BODY
#  Parse the command line switches
    optioncount = len(sys.argv)
    SAC = False
    outputfile_defined = False
    filelist = []
    dir=""
    extension = '.txt'
    if optioncount > 1:

        if optioncount == 4:
            if sys.argv[3] == '-s':
                SAC = True
#               print "SAC set to true."
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

        
        for n in range(len(filelist)):
            if extension in filelist[n]:
                if len(filelist)>1:
                    infile = dir+"/"+filelist[n]

                if string.find(infile,'.') > 0:
                    outfile = infile[:string.find(infile,'.')]+'.sac'
                else:
                    outfile = infile +'.sac'

                PNE = load(infile)

                Comment = PNE[0][0][1]
                Stname = PNE[0][1][1][:7]
                Component = PNE[0][2][1][:3]
                St_time = time.strptime(PNE[0][4][1][:-4],"%d_%b_%Y_%H:%M:%S")
                Frac_sec = int(PNE[0][4][1][21:])
                CF = np.float32(PNE[0][5][1])

                TC = PNE[0][6][1]
                Offset = float(PNE[1][0][0])
#
#                       Delta is calculated from the offset time of last sample 
#                       minus offset of first sample / total number of samples.

            
                Delta = (float(PNE[1][len(PNE[1])-1][0])-Offset)/(len(PNE[1])-1)

#                       Load Data array
#                       Samples in file are multiplied by 10,000 to convert from
#                       measurements of centimeters to microns, then it's divided by
#                       the Amplification (conversion) factor, known as CF

                Data = []
                for n in range (len(PNE[1])-1):
                    Datum = np.float32(np.float32(PNE[1][n][1])*10000.0/CF)
                    Data.append(Datum)
                t = SacIO()
                b = np.arange(len(Data),dtype=np.float32)
                for n in range(len(Data)): #   Load the array with time-history data
                    b[n] = Data[n]
                t.fromarray(b)

                t.SetHvalue('delta', Delta)
                t.SetHvalue('nzyear',St_time.tm_year)
                t.SetHvalue('nzjday',St_time.tm_yday)
                t.SetHvalue('nzhour',St_time.tm_hour)
                t.SetHvalue('nzmin',St_time.tm_min)
                t.SetHvalue('nzsec', St_time.tm_sec)
                t.SetHvalue('nzmsec', Frac_sec)
                t.SetHvalue('kstnm',Stname[:7])
                t.SetHvalue('kcmpnm',Component)
                t.WriteSacBinary(outfile)

                print "File written to {}".format(outfile)
            else:
                print "{} not processed.".format(filelist[n])
        
    else:
        print "Useage: PNE2SAC infile.txt (outfile.asc)"
        print "Or, PNE2SAC target_directory target_extension(like .txt)"
        print "No infile or directory specified."
        print len(sys.argv)

#
# Check and run the main function here:
#
if __name__ == '__main__':
  main()
