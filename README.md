Converters
==========

Seismic Datafile converters to SAC or miniseed
The converters found within this directory are used for managing seismic data
from various local networks throughout the Russian federation. These networks use many different types of seismometers
as well as digitizing systems, all of which employ unique data formats. One of the more common systems we run across
is manufactured by Symmetric Research, and it uses two different data formats: A DAT format (various flavors), and an
OUT format, which is an adaptation created by the Magadan / Yakutia regional networks. (DAT2SAC, OUT2SAC) 

Other formats we deal with include a hand-digitizng waveform tracing tool that is used for transcribing analog records
into a digital format for the preservation of historical seismograms. This is PNE2SAC. It will take any ascii file that
conforms to the data format within, and turn it into a SAC file. 

Lastly, there are also converters that convert certain comma separated variable text file formats into SAC files.

The files, once converted into SAC, should be readable by any ObsPy application, and likely, any SAC utilities. 
