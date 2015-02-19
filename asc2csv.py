__author__ = "Daniel Burk <burkdani@msu.edu>"
__version__ = "20140814"
__license__ = "MIT"



import os
import sys
                                           #from scipy import signal
                                           #from scipy.integrate import simps
                                           #import pylab as plt
                                           #import numpy as np
                                           #import scipy as sp
import csv
import calendar
import string                              # We use the string.find utility for file naming
import time                                # we use the sleep function to enable conversion of the DAT
# import subprocess                          # used for the execution of command-line programs
                                           # from obspy.core import read, Trace, Stream, UTCDateTime
                                           # from obspy.sac import SacIO


class asc2csv(object):
    '''asc2csv.py is a utility for batch converting Symmetric Research ASC files into csv files
       for a whole directory. It has only one command line
       switch, and that is to point it at a directory for converting the dat files.

       The program will convert asc files into separate csv files that match the format of the
       DAT2ASC.EXE program from symmetric research.
        
       

       Syntax: asc2csv target_directory  

       

       Typical useage:
       <ObsPy> C:\Python27\scripts> python asc2csv.py c:/calibration/station/ 

    '''



def convert(infile):
    print 'File name is {}'.format(infile)
    with open(infile,'r') as fin:
        list = csv.reader(fin,delimiter = ' ')
        rowcnt = 0
    
        stack = []
        header = []
        for row in list:
        
            if len(row)== 6:        # all headers are placed within the header list
                header.append(row)
    
        
            else:
                if row:
                    stack.append(row) # all data is placed in the stack
        
                                # Parse the data into a useful list:
                                # first element: Channel 0 data
                                # second element: Channel 1 data
                                # third element: Channel 2 data
                                # and so on to nth data
                                # followed by the time of day for each sample.
    
    Sttime = header[0][0]                          # Retrieve start time text field
    interval = str(float(header[0][4])/1000)       # Retrieve the sample interval
    numsamples = int(header[0][5])
    t = time.strptime(Sttime[0:14],"%Y%m%d%H%M%S") # convert the first time into linux time
    starttime = calendar.timegm(t)                 # convert the time into linux time
    starttime_r = float(Sttime[14:17])/1000        # Get the remainder of fractions of second
    sampletime = []
    sampletime_r = []
    for n in range(0,numsamples):      # Create sampletime and Fracsec lists
        ET = (starttime_r+(n * float(interval)))
        sampletime_r.append('%.3f' % (ET-int(ET)))
        sampletime.append(starttime+int(ET))

    data = [[],[],[],[]]                              # currently handles only four channels

    for n in range(0,len(header)):                    # Break the stack into it's requisite signal channels
        for i in range(0,numsamples):
            data[n].append(stack[int(n*numsamples)+i])

                                                  # Next, convert the file to the existing csv format
                                                  # First, populate the header with the appropriate information
    head = []
    for i in range(0,16):
        head.append("na")
    head[0] = "Pt#"
    head[1] = header[0][2]          # put the channel names into the head
    head[2] = header[2][2]
    head[3] = header[2][2]
    head[4] = header[1][2]
    
    head[12]="Linuxtime"            # Add the descriptor for what the time column represents
    head[15]="FracSec"              # last column is the fraction of the second paired with linux time


    
    outfile = infile[string.rfind(infile,"\\")+1:string.find(infile,'.')]+".csv"

    row = ['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']
    with open(outfile,'wb') as csvfile: # use 'wb' in place of 'a' if you want to overwrite the file.
                outrow = csv.writer(csvfile, delimiter = ",",
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
                outrow.writerow(head)
                for i in range(0,numsamples):
                    row[1] = int(float(data[0][i][0]))   # x channel (function generator)
                    row[2] = int(float(data[2][i][0]))   # z channel (signal coil)
                    row[3] = int(float(data[2][i][0]))   # z channel (signal coil)
                    row[4] = int(float(data[1][i][0]))   # y channel (laser position sensor) 
                    row[12] = sampletime[i]
                    row[15] = sampletime_r[i]
                    outrow.writerow(row)





def main():
                                      #           MAIN PROGRAM BODY
                                      #  Parse the command line switches
                                      # Commmand example: c:\Python27>Python.exe Sigcal.py c:\seismo\caldata\momo
                                      # where momo is the working directory containing the csv files
                                      # as well as the calibration control file, c:\seismo\caldta\calcontrol.csv
                                      # The third option can designate an optional location for the calcontrol file.
                                      #
    optioncount = len(sys.argv)
    outputfile_defined = False
    filelist = []
    dir=""
    infile = ""
    directory = os.getcwd()
    directory = directory.replace("/","\\")
    print "This is the directory name as taken from the computer:",directory
    if optioncount > 1:
        
        
        directory = sys.argv[1]
        directory = directory.replace("/","\\")
        if directory[-1:] !="\\":
                directory = directory+"\\"
        try:
            filelist = os.listdir(directory)            
        except:
            print "Command line parameter must consist of a valid directory. {}".format(directory)
            sys.exit(0)
    else:
        filelist = os.listdir(directory)
    

    for n in range(0,len(filelist)):                                
        if ".asc" in filelist[n]:
            
            convert(filelist[n])

#
# Check and run the main function here:
#
if __name__ == '__main__':
  main()
 