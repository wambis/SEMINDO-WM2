#!/usr/bin/env python
# read_csem.py

"""
Packages for reading eod/h files
Developed by Haydar Karaoglu (May, 2015)
Modified by Federico D. Munch (April, 2020)
""" 

# Load python packages
import os
import numpy as np
import struct
import re

# ---------------------------------------------------------------------- 
# EODH Class - eismic data
# ----------------------------------------------------------------------     
class eohPacket:    

    def __init__(self,f=None):
        """
        Initialize seismic data header -- if f is specified it must
        be a file object and is used to fill the structure, 
        otherwise we get a blank structure which we can fill
        """
        
        #self._fmt = "1i12c4c2f2i2f4c" # Old Haydar formatting
        self._fmt = "1i4s4s4s1s1s1s1s2f2i2f4s"
        self._size = struct.calcsize(self._fmt)
       
        if f != None:
            self.read(f)
            
    def read(self,f):
        """
        read wavepacket header from file object
        """
        
        # record file and point at which this wavepacket header is stored
        self.eoh_name = os.path.basename(f.name)
        self.eoh_loc = f.tell()
        self.event = os.path.basename(f.name).replace(".eoh","")
        self.eod_name = os.path.basename(f.name).replace(".eoh",".eod")
        
        read_data = f.read(self._size)
        if (len(read_data) != self._size):
            raise IOError("Struct_packhdrbr_st: not enough data to fill class")

        packhdrbr_st = struct.unpack(self._fmt,read_data)
        self.rcno = int(packhdrbr_st[0])
        self.netwk = "".join(item[0] for item in packhdrbr_st[1:5])
        self.chnl  = "".join(item[0] for item in packhdrbr_st[5:9])
        self.stn   = "".join(item[0] for item in packhdrbr_st[9:13])
        self.compnt = packhdrbr_st[13]
        self.dttype = packhdrbr_st[14]
        self.reftime = packhdrbr_st[15]
        self.padd = packhdrbr_st[16]
        self.starttm = float(packhdrbr_st[17])
        self.smplintv = float(packhdrbr_st[18])
        self.ndata = int(packhdrbr_st[19])
        self.locatn = int(packhdrbr_st[20])
        self.slat   = float(packhdrbr_st[21])
        self.slon   = float(packhdrbr_st[22])
        self.extra  = "".join(item[0] for item in packhdrbr_st[23:27])
        
        """
        Return name of CSEM ascii file given station and component name
        """
        stn = re.sub('[( ]', '', self.stn)  #Erase empty space in station name (Identified)
        netwk = re.sub('[( ]', '', self.netwk)  #Erase empty space in network name (Maynot exist)
                
        if len(stn) > 4:
            raise Exception("eoh: station name length error")
        else:
            n = 4 - len(stn)
            stn = "_"*n + stn
            
        if len(netwk) > 4:
            raise Exception("eoh: network name length error")
        else:
            n = 4 - len(netwk)
            stn = "_"*n + netwk + '.' + stn     

        if self.compnt not in ["Z","T","L"]:
            raise Exception("eoh: component name error")

        self.longstn = "U{0}_{1}".format(self.compnt,stn)
        self.csem = False
        
        
#        self.rcno = 0              #int rcno
#        self.netwk  = 'aaaa'       #char netwk[4]
#        self.chnl   = 'aaaa'       #char chnnl[4];
#        self.stn    = 'aaaa'       #char stn[4];
#        self.compnt = 'F'          #char compnt;        /* N=north; T=transvers etc. */
#        self.dttype = 'a'          #char dttype;        /* a=accel., c=counts  etc. */
#        self.reftime= 'c'          #char reftime;       /* c=cmt time; b=bulletin time*/
#        self.padd   = '0'          #char padd;
#        self.starttm=  0.0         #float starttm;      /* in seconds from reftime */
#        self.smplintv= 0.0         #float smplintv;     /* in seconds */
#        self.ndata   = 0;          #int ndata;
#        self.locatn  = 0;          #int locatn;        /* starting locatn in *.eod in
#                                   #                    * bytes; first one is 0 
#                                   #                    */
#        self.slat = 0.0;           #float slat;          /*station lat in degrees */
#        self.slon = 0.0;           #float slon;          /*station lon in degrees */
#        self.extra= 'aaaa';        #char extra[4];
        
    def write(self,f,offset=None):
        """
        write eohd header object to file (object)
        if offset is specified then we first seek to the 
        specified offset from the start of the file
        """
        
        if offset != None:
            f.seek(offset,0)
        
        header = struct.pack(self._fmt,self.rcno,self.netwk,self.chnl,\
        self.stn,self.compnt,self.dttype,self.reftime,self.padd,\
        self.starttm,self.smplintv,self.ndata,self.locatn,self.slat,\
        self.slon,self.extra)

        f.write(header)

class eodpacket:
    def __init__(self):
        self.header = None
        self.data = None

    def read(self,f,ndata):
        self.header = True
        self.data = np.fromfile(f,dtype=np.float32,count=ndata)

        
class eodData:

    def __init__(self,eohfile,datafile=None):
        self._header_file = eohfile
        
        # if datafile is specified, then we set the self._data_file to that
        # otherwise the data file is expected to be in the same directory
        # as the wph file
        
        if datafile != None:
            self._data_file = datafile
        else:
            self._data_file = re.sub('eoh$','eod',self._header_file)

        
        self.event = None
        self.header_read = False
        self.data_read = False

    def __check_header_loaded(self):
        if self.header_read == False:
            raise Exception("eoData instance has not had header loaded")

    def __check_data_loaded(self):
        if self.data_read == False:
            raise Exception("eoData instance had not had data loaded")

    def read_header(self):
        if not os.path.isfile(self._header_file):
            raise IOError("File %s does not exist" % self._header_file)
       
        # open eoh file for binary read
        eoh = open(self._header_file,'rb')
        self.offset = eoh.read(16) # Offset at the beginning of header files

        # read event details
        self.event = self._header_file[-12:-4]

        # read trace details - Struct_tracehdr_st class raises exception
        # when end of file is reached
        self.eod_hdrs = []
        try:
            while (1):
                eodhdr_st = eohPacket(eoh)
                self.eod_hdrs.append(eodhdr_st)
        except IOError:
            pass
        
        self.neod = len(self.eod_hdrs)
        eoh.close()
        self.header_read = True
        
    def write_header(self,filename):
        
        f=open(filename,'wb')        
        f.write(self.offset.encode())

        for eohid in range(len(self.eod_hdrs)):
            header=self.eod_hdrs[eohid]    
            
            if header.csem == True:   # If it is not csem synthetics, do not write to the new file             

                header.netwk = str(header.netwk).encode()
                header.chnl = str(header.chnl).encode()
                header.stn = str(header.stn).encode()
                header.compnt = str(header.compnt).encode()
                # Old Haydar formatting
                """
                data = struct.pack(header._fmt,header.rcno,
                                   header.netwk[0],header.netwk[1],header.netwk[2],header.netwk[3],
                                   header.chnl[0],header.chnl[1],header.chnl[2],header.chnl[3],
                                   header.stn[0],header.stn[1],header.stn[2],header.stn[3],
                                   header.compnt,
                                   header.dttype.encode(),
                                   header.reftime.encode(),
                                   header.padd.encode(),
                                   header.starttm,
                                   header.smplintv,header.ndata,header.locatn,header.slat,header.slon,
                                   header.extra[0].encode(),
                                   header.extra[1].encode(),
                                   header.extra[2].encode(),
                                   header.extra[3].encode() )
                """
                data = struct.pack(header._fmt,header.rcno,
                                   header.netwk,
                                   header.chnl,
                                   header.stn,
                                   header.compnt,
                                   header.dttype.encode(),
                                   header.reftime.encode(),
                                   header.padd.encode(),
                                   header.starttm,
                                   header.smplintv,header.ndata,header.locatn,header.slat,header.slon,
                                   str(header.extra).encode() )

                f.write(data)
            
        f.close()

    def read_data(self):
        # read header if not done separately
        if self.event == None:
            self.read_header()

        # check data file exists
        if not os.path.isfile(self._data_file):
            raise IOError("File %s does not exist" % self._data_file)

        eod = open(self._data_file,'rb')
        self.eod_data = []

        for eodid in self.eod_hdrs:
            eod.seek(eodid.locatn)
            try:
                    ed = eodpacket()
                    ed.read(eod,eodid.ndata)
                    self.eod_data.append(ed)
            except IOError:
                pass
            
        self.data_read = True
        
    def write_data(self,filename):
        
        f=open(filename,'wb')        

        for eohid in range(len(self.eod_hdrs)):
            header=self.eod_hdrs[eohid]    
            if header.csem == True:   # If it is not csem synthetics, do not write to the new file  
                data=self.eod_data[eohid].data
                outdata = np.array(data,dtype=np.float32)
                outdata.tofile(f)
            
        f.close()

def get_eodh_file_lists(eodDataDir):
    """given the directory containing eod/eph files return lists containing
    their paths"""

    eodDirList = [ f for f in os.listdir(eodDataDir) 
                  if os.path.isfile(os.path.join(eodDataDir,f)) ]
    eohFiles = [ f for f in eodDirList if re.match(".*\.eoh$",f) ]
    eohFiles.sort()

    eohMatches = []
    eodMatches = []

    for eohFile in eohFiles:
        eohPath = os.path.join(eodDataDir,eohFile)
        eodFile = eohFile.replace(".eoh",".eod")
        eodPath = os.path.join(eodDataDir,eodFile)

        if os.path.isfile(eodPath):
            eohMatches.append(eohPath)
            eodMatches.append(eodPath)
      
    return eohMatches, eodMatches

       
        
        
        
