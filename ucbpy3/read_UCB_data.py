#!/usr/bin/env python

import math
import numpy as np
import os
import re
import struct
import sys

# ********************************************************************** #
# Classes representing C structs
# ********************************************************************** #
class Struct_packhdrbr_st:
    def __init__(self,f=None):
        """
        Initialize wavepacket header -- if f is specified it must
        be a file object and is used to fill the structure, 
        otherwise we get a blank structure which we can fill
        """
        self._fmt = "=3h2s9f"
        self._size = struct.calcsize(self._fmt)
   
        # create empty object
        self.wpd_name = None
        self.wpd_loc = None
        self.phase = 0
        self.ndata = 0  
        self.id = 0
        self.t0 = 0. 
        self.gv1 = 0.
        self.gv2 = 0. 
        self.pv1 = 0.
        self.pv2 = 0.
        self.rmsd = 0.
        self.rmsr = 0.
        self.rmss = 0.
        self.weight = 0.

        # read in if we've specified a file object
        if f != None:
            self.read(f)

    def read(self,f):
        """
        read wavepacket header from file object
        """
        
        # record file and point at which this wavepacket header is stored
        self.wpd_name = os.path.basename(f.name)
        self.wpd_loc = f.tell()

        read_data = f.read(self._size)
        if (len(read_data) != self._size):
            raise IOError("Struct_packhdrbr_st: not enough data to fill class")

        packhdrbr_st = struct.unpack(self._fmt,read_data)
        self.phase = int(packhdrbr_st[0])
        self.ndata = int(packhdrbr_st[1])
        self.id = int(packhdrbr_st[2])
        self.t0 = float(packhdrbr_st[4])
        self.gv1 = float(packhdrbr_st[5])
        self.gv2 = float(packhdrbr_st[6])
        self.pv1 = float(packhdrbr_st[7])
        self.pv2 = float(packhdrbr_st[8])
        self.rmsd = float(packhdrbr_st[9])
        self.rmsr = float(packhdrbr_st[10])
        self.rmss = float(packhdrbr_st[11])
        self.weight = float(packhdrbr_st[12])
        
    def write(self,f,offset=None):
        """
        write wavepacket header object to file (object)
        if offset is specified then we first seek to the 
        specified offset from the start of the file
        """

        if offset != None:
            f.seek(offset,0)
       
        dummy = "00"  # there is padding in the structure
        data = struct.pack(self._fmt,self.phase,self.ndata,self.id,dummy,
                           self.t0,self.gv1,self.gv2,self.pv1,self.pv2,
                           self.rmsd,self.rmsr,self.rmss,self.weight)

        f.write(data)

    def write_header_and_data(self,f,data,offset=None,update=True):
        """
        Write wavepacket header object to file (object), followed by 
        data (will check consistency with header).   
        Updates rmsd in header to rms value of new data.
        If offset is specified we first seek to the specified offset
        from the start of the file (before writing header)
        """

        if data.size != self.ndata:
            print("Data array size: %d,  packhdr ndata: %d" %(data.size,self.ndata))
            print("Data t0: %f" % self.t0)
            print("Trying to write new data: ")
            print (data)
            raise IOError("Struct_packhdr_st.write_header_and_data - data array length inconsistent with ndata")

        # update rmsd based on new data
        if update:
            self.rmsd = np.sqrt(np.mean(data**2)) 

        # write header then data
        self.write(f,offset)
        outdata = np.array(data,dtype=np.float32)
        outdata.tofile(f)
#        f.write(outdata.tobytes())


class Struct_tensor_st:
    def __init__(self,rr,tt,pp,rt,rp,tp):
        self.rr = float(rr)
        self.tt = float(tt)
        self.pp = float(pp)
        self.rt = float(rt)
        self.rp = float(rp)
        self.tp = float(tp)

class Struct_title_st:
    def __init__(self,f):
        self._fmt = "=8s10f"
        self._size = struct.calcsize(self._fmt)
        
        read_data = f.read(self._size)
        if (len(read_data) != self._size):
            raise IOError("Struct_title_st: not enough data to fill class")

        title_st = struct.unpack(self._fmt,read_data)
        self.event = title_st[0]
        self.theta = float(title_st[1])
        self.phi = float(title_st[2])
        self.depth = float(title_st[3])
        self.moment = Struct_tensor_st(title_st[4],title_st[5],title_st[6],
                                       title_st[7],title_st[8],title_st[9])
        self.dt = float(title_st[10])

class Struct_tracehdr_st:
    def __init__(self,f):
        self._fmt = "=4si10fh2c4s4s4s"
        self._size = struct.calcsize(self._fmt)
        
        read_data = f.read(self._size)
        if (len(read_data) != self._size):
            raise IOError("Struct_tracehdr_st: not enough data to fill class")

        # n.b. for chnnl and extr, we get rid of any non-alphanumeric characters
        # before assignment, as the C struct doesn't use null terminated strings
        # so if we don't do this we get garbage...
        tracehdr_st = struct.unpack(self._fmt,read_data)
        self.stn = tracehdr_st[0].strip()
        self.locatn = int(tracehdr_st[1])
        self.theta = float(tracehdr_st[2])
        self.phi = float(tracehdr_st[3])
        self.delta = float(tracehdr_st[4])
        self.az = float(tracehdr_st[5])
        self.dip = float(tracehdr_st[6])
        self.smplintv = float(tracehdr_st[7])
        self.w1 = float(tracehdr_st[8])
        self.w2 = float(tracehdr_st[9])
        self.w3 = float(tracehdr_st[10])
        self.w4 = float(tracehdr_st[11])
        self.id = int(tracehdr_st[12])
        self.reftime = tracehdr_st[13]
        self.comp = tracehdr_st[14]
        self.netwk = tracehdr_st[15].strip()
        self.chnnl = re.sub(r'\W','',tracehdr_st[16].strip())
        self.extr = re.sub(r'\W','',tracehdr_st[17].strip())
        
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

        if self.comp not in ["Z","T","L"]:
            raise Exception("wph: component name error")

        self.longstn = "U{0}_{1}".format(self.comp,stn)


class Struct_tracehdrH_st(Struct_tracehdr_st):
    """ this subclasses Struct_tracehdr_st since there is only 
    one additional variable, which is written at the end of the 
    struct """
    def __init__(self,f):
        Struct_tracehdr_st.__init__(self,f)
        self.locatnA = struct.unpack('=i',f.read(struct.calcsize('=f'))) 

# ********************************************************************** #
# classes to read data files
# ********************************************************************** #
class Wavepacket:
    def __init__(self):
        self.header = None
        self.data = None

    def read(self,f):
        self.header = Struct_packhdrbr_st(f)
        self.data = np.fromfile(f,dtype=np.float32,count=self.header.ndata)


class wpData:

    def __init__(self,wphfile,sem=False,datafile=None):
        self._header_file = wphfile
        
        # if datafile is specified, then we set the self._data_file to that
        # otherwise the data file is expected to be in the same directory
        # as the wph file, and sem==True wil expect it to have .wph replaced
        # by .smd, sem==False with .wpd
        if datafile != None:
            self._data_file = datafile
        elif sem == False:
            self._data_file = re.sub('wph$','wpd',self._header_file)
        else:
            self._data_file = re.sub('wph$','smd',self._header_file)
        
        self.event = None
        self.trace_hdrs = []
        self.traces = []
        self.trace_nwp = []
        self.ntrace = 0
        self.header_read = False
        self.data_read = False

    def __check_header_loaded(self):
        if self.header_read == False:
            raise Exception("wpData instance has not had header loaded")

    def __check_data_loaded(self):
        if self.data_read == False:
            raise Exception("wpData instance had not had data loaded")

    def read_header(self):
        if not os.path.isfile(self._header_file):
            raise IOError("File %s does not exist" % self._header_file)
       
        # open wph file for binary read
        wph = open(self._header_file,'rb')

        # read event details
        self.event = Struct_title_st(wph)

        # read trace details - Struct_tracehdr_st class raises exception
        # when end of file is reached
        self.trace_hdrs = []
        try:
            while (1):
                tracehdr_st = Struct_tracehdr_st(wph)
                self.trace_hdrs.append(tracehdr_st)
        except IOError:
            pass
        self.ntrace = len(self.trace_hdrs)
        wph.close()
        self.header_read = True

    def read_data(self):
        # read header if not done separately
        if self.event == None:
            self.read_header()

        # check data file exists
        if not os.path.isfile(self._data_file):
            raise IOError("File %s does not exist" % self._data_file)

        self.traces = []
        self.trace_nwp = []
        wpd = open(self._data_file,'rb')
        for thdr in self.trace_hdrs:
            wpd.seek(thdr.locatn)
            wpackets = []
            try:
                while(1):
                    wp = Wavepacket()
                    wp.read(wpd)
                    if wp.header.id != thdr.id:
                        break
                    else:
                        wpackets.append(wp)
            except IOError:
                pass
            self.traces.append(wpackets)
            self.trace_nwp.append(len(wpackets))
        self.data_read = True

    def get_event(self):
        return self.event.event

    def get_event_coords(self):
        sc = zip([self.event.event], [90. - math.degrees(self.event.theta)], [math.degrees(self.event.phi)],[self.event.depth])
        return sc

    def get_stns_comp_coords(self):
        self.__check_header_loaded()
        stations = []
        comps = []
        lats = []
        lons = []
        for thdr in self.trace_hdrs:
            lat = 90. - math.degrees(thdr.theta)
            lon = math.degrees(thdr.phi)
            stations.append(thdr.stn)
            comps.append(thdr.comp)
            lats.append(lat)
        lons.append(lon)
        sc = zip(stations,comps,lats,lons)
        return sc	

    def get_stns(self):
        self.__check_header_loaded()
        stations = []
        for thdr in self.trace_hdrs:
            stations.append(thdr.stn)
        stations.sort()
        return stations

    def get_stns_comp(self):
        self.__check_header_loaded()
        stations = []
        comps = []
        for thdr in self.trace_hdrs:
            stations.append(thdr.stn)
            comps.append(thdr.comp)

        sc = zip(stations,comps)
        return sc

    def get_nets_stns(self,sort='net'):
        self.__check_header_loaded()
        stations = []
        networks = []
        for thdr in self.trace_hdrs:
            stations.append(thdr.stn)
            networks.append(thdr.netwk)

        # sort returned list, either by network then station ('net')
        # or by station then network ('stn')
        ns = zip(networks,stations)
        if sort == 'net':
            cols = (0,1)
        elif sort == 'stn':
            cols = (1,0)
        else:
            raise ValueError("Incorrect value of parameter sort")

        return sort_zipped_list(ns,cols)

    def get_deltas(self,sort='delta',radians=False):
        stations = []
        delta = []
        for thdr in self.trace_hdrs:
            stations.append(thdr.stn)
            if radians == True:
                delta.append(thdr.delta)
            else:
                delta.append(math.degrees(thdr.delta))

        # return sorted zipped lists
        sd = zip(stations,delta)
        if sort == 'delta':
            cols = (1,0)
        elif sort == 'stn':
            cols = (0,1)
        else:
            raise ValueError("Incorrect value of parameter sort")

        return sort_zipped_list(sd,cols)

    def print_metadata(self,sort='stn'):
        self.__check_header_loaded()
        self.__check_data_loaded()

        # Determine sort order
        index = []
        sort_var_1 = []
        sort_var_2 = []

        for i,thdr in enumerate(self.trace_hdrs):
            index.append(i)
            if sort == 'stn':
                sort_var_1.append(thdr.stn)
                sort_var_2.append(thdr.netwk)
            elif sort == 'net':
                sort_var_1.append(thdr.netwk)
                sort_var_2.append(thdr.stn)
            elif sort == 'delta':
                sort_var_1.append(thdr.delta)
                sort_var_2.append(thdr.stn)

        slist = zip(index,sort_var_1,sort_var_2)
        slist = sort_zipped_list(slist,(1,2))
        ind_order = zip(*slist)[0]

        for i in ind_order:
            thdr = self.trace_hdrs[i]
            net = thdr.netwk.strip().replace('(','')
            stn = thdr.stn.strip()
            delta = math.degrees(thdr.delta)
            lat = 90. - math.degrees(thdr.theta)
            lon = math.degrees(thdr.phi)

            nwp = self.trace_nwp[i]
            phases = []
            for j in range(nwp):
                phases.append(str(self.traces[i][j].header.phase))

            print("%-4s %-4s %8.4f %8.4f %9.4f %2d %s" % \
                    (net, stn, delta, lat, lon, nwp, ' '.join(phases)))


# ********************************************************************** #
# utility functions
# ********************************************************************** #
def sort_zipped_list(zlist,cols):
    """
    sort zipped list zlist by columns in the order specified in the 
    tuple cols -- i.e. (1,0) sorts by column 1 then column 0
    """

    for col in reversed(cols):
        zlist.sort(key=lambda x: x[col])
    return zlist


