#!/usr/bin/env python
# **********************************************************************
# data_selector.py
# **********************************************************************
# Routines to read wavepacket data into a single structure for
# each wavepacket, combining real and SEM data, and then to 
# select subsets of the data.   
#
# TODO: probably combine this into the read_UCB_data module
#
# JB, Sept 2012
# **********************************************************************

import pickle as cPickle
import math
import os
import re
import sys
import warnings

import numpy as np
import scipy.interpolate
import scipy.stats
import scipy.signal
import read_UCB_data
import obspy.signal.polarization as pol
import obspy.signal.tf_misfit as misfit

# ---------------------------------------------------------------------- 
# WavpacketData Class - combines all wavepacket info (data and sem)
# ---------------------------------------------------------------------- 
class WavepacketData:
    """Class combines all data for a wavepacket into a single object,
    including both real and sem data"""

    def __init__(self):
        self._header_populated = False
        self._data_populated = False
        self._data_env_populated = False
        self._sem_populated = False
        self._sem_env_populated = False

    def import_data(self,event_hdr,trace_hdr,wpkt):
        if self._header_populated == True:
            warnings.warn("Wavepacket already loaded - new data will replace it!")

        # event parameters from event header
        self.event_id = event_hdr.event
        
        self.event_lat = 90. - math.degrees(event_hdr.theta)
        self.event_lon = math.degrees(event_hdr.phi)
        self.event_depth = event_hdr.depth
        self.event_dt = event_hdr.dt

        self.moment_rr = event_hdr.moment.rr
        self.moment_tt = event_hdr.moment.tt
        self.moment_pp = event_hdr.moment.pp
        self.moment_rt = event_hdr.moment.rt
        self.moment_rp = event_hdr.moment.rp
        self.moment_tp = event_hdr.moment.tp

        # trace parameters from trace header
        self.stn = trace_hdr.stn
        self.stn_lat = 90. - math.degrees(trace_hdr.theta)
        self.stn_lon = math.degrees(trace_hdr.phi)
        self.delta = trace_hdr.delta
        self.az = trace_hdr.az
        self.dip = trace_hdr.dip
        self.smplintv = trace_hdr.smplintv
        self.filt_w1 = trace_hdr.w1
        self.filt_w2 = trace_hdr.w2
        self.filt_w3 = trace_hdr.w3
        self.filt_w4 = trace_hdr.w4
        self.reftime = trace_hdr.reftime #n.b. just a single char, not actual ref time
        self.comp = trace_hdr.comp
        self.netwk = trace_hdr.netwk
        self.chnnl = trace_hdr.chnnl
        self.extr = trace_hdr.extr
        self.longstn = trace_hdr.netwk+'.'+trace_hdr.stn

        # wavepacket parameters from wavepacket header
        self.wpd_name = wpkt.header.wpd_name
        self.wpd_loc = wpkt.header.wpd_loc
        self.phase = wpkt.header.phase
        self.ndata = wpkt.header.ndata
        self.t0 = wpkt.header.t0
        self.gv1 = wpkt.header.gv1
        self.gv2 = wpkt.header.gv2
        self.pv1 = wpkt.header.pv1
        self.pv2 = wpkt.header.pv2
        self.rmsd_data = wpkt.header.rmsd
        self.rmsr_data = wpkt.header.rmsr
        self.rmss_data = wpkt.header.rmss
        self.weight_data = wpkt.header.weight
        
        self.u_rmss = self.rmss_data
        self.u_rmsd = self.rmsd_data
        self.u_rmsr = self.rmsr_data

        self.env_rmss = 0.
        self.env_rmsd = 0.
        self.env_rmsr = 0.
        
        self.sp_amp_rmsr = 0.
        self.sp_amp_rmsd = 0.
        self.sp_amp_rmss = 0.
        
        self.sp_phase_rmsr = 0.
        self.sp_phase_rmsd = 0.
        self.sp_phase_rmss = 0.
        
        self.sp_amp_correlation = 0.
        self.sp_phase_correlation = 0.
        
        self.correlation = 0.
        self.iphase_correlation = 0.
        self.spec_correlation = 0.
        self.env_correlation = 0.
        self.l2_misfit = 0.
        self.l2_misfit_normalized = 0.
        self.l1_misfit = 0.
        self.l1_misfit_normalized = 0.
        self.insfreq_correlation = 0.
        self.pgf = -1
        self.egf = -1
        self.globalmisfit = 0
        self.env_globalmisfit = 0
        
        self.data_smoothness_factor = 0
        self.sem_smoothness_factor = 0

        
        # wavepacket data
        self.data = wpkt.data
        self.__calculate_data_env()

        # consistency checks - store id for checking against sem data when
        # read in, but we don't otherwise need it so it's an underscore variable
        if trace_hdr.id != wpkt.header.id:
            raise IOError("Internal inconsistency in data read: id")
        else:
            self._id = trace_hdr.id

        # time array to go with data array
        self.t = np.arange(self.ndata)*self.smplintv

        # data structure now populated
        self._header_populated = True
        self._data_populated = True

    def import_sem(self,wpkt):
        """ add SEM data to structure, and check its wavepacket header
        for consistency """

        if self._data_populated == False:
            raise RuntimeError("SEM data may not be imported before real data")
        if self._sem_populated == True:
            warnings.warn("SEM data already loaded - new data will replace it!")

        # check for consistency
        if wpkt.header.id != self._id:
            raise IOError("Internal inconsistency in SEM data read: id")
        if wpkt.header.phase != self.phase:
            raise IOError("Internal inconsistency in SEM data read: phase")
        if wpkt.header.ndata != self.ndata:
            raise IOError("Internal inconsistency in SEM data read: ndata")
        if wpkt.header.t0 != self.t0:
            raise IOError("Internal inconsistency in SEM data read: t0")
        if wpkt.header.gv1 != self.gv1:
            raise IOError("Internal inconsistency in SEM data read: gv1")
        if wpkt.header.gv2 != self.gv2:
            raise IOError("Internal inconsistency in SEM data read: gv2")
        if wpkt.header.pv1 != self.pv1:
            raise IOError("Internal inconsistency in SEM data read: pv1")
        if wpkt.header.pv2 != self.pv2:
            raise IOError("Internal inconsistency in SEM data read: pv2")
       
        # assign additional header variables
        self.smd_name = wpkt.header.wpd_name
        self.rmsd_sem = wpkt.header.rmsd
        self.rmsr_sem = wpkt.header.rmsr
        self.rmss_sem = wpkt.header.rmss
        self.weight_sem = wpkt.header.weight

        # assign data
        self.sem_data = wpkt.data
        self.__calculate_sem_env()
        self._sem_populated = True

        # calculate various correlations, misfits, etc between sem and data
        # (including between envelopes)
        self.__calculate_correlation_coeff()
#        self.__calculate_spec_correlation_coeff()
#        self.__calculate_iphase_correlation_coeff()
#        self.__calculate_insfreq_correlation_coeff()
#        self.__calculate_cross_correlation()
        self.__calculate_l1_misfit()
        self.__calculate_l2_misfit()
        self.__calculate_env_correlation_coeff()
#        self.__calculate_env_cross_correlation()
#        self.__calculate_env_l1_misfit()
        self.__calculate_env_l2_misfit()
        self.__calculate_global_misfit()
        
        self.__calculate_spec_correlation_coeff()

#        self.__calculate_gof()

    def return_packhdr(self):
        """
        return wavepacket header instance with values from this
        class.  This can be used if changes (e.g. weighting) are 
        made which need to be written to the wpd files (as the 
        Struct_packhdrbr_st class now has a write method which
        can write to a specified point in the wpd file).
        """

        if self._header_populated != True:
            return None

        pkhdr = read_UCB_data.Struct_packhdrbr_st()
        pkhdr.wpd_name = self.wpd_name
        pkhdr.wpd_loc = self.wpd_loc
        pkhdr.phase = self.phase
        pkhdr.ndata = self.ndata
        pkhdr.id = self._id
        pkhdr.t0 = self.t0
        pkhdr.gv1 = self.gv1
        pkhdr.gv2 = self.gv2
        pkhdr.pv1 = self.pv1
        pkhdr.pv2 = self.pv2
        pkhdr.rmsd = self.rmsd_data
        pkhdr.rmsr = self.rmsr_data
        pkhdr.rmss = self.rmss_data
        pkhdr.weight = self.weight_data

        return pkhdr

    def return_packhdr_sem(self):
        """
        return wavepacket header instance with values from this
        class.  This can be used if changes (e.g. weighting) are 
        made which need to be written to the wpd files (as the 
        Struct_packhdrbr_st class now has a write method which
        can write to a specified point in the wpd file).
        """

        if self._header_populated != True:
            return None

        pkhdr = read_UCB_data.Struct_packhdrbr_st()
        pkhdr.wpd_name = self.smd_name
        pkhdr.wpd_loc = self.wpd_loc
        pkhdr.phase = self.phase
        pkhdr.ndata = self.ndata
        pkhdr.id = self._id
        pkhdr.t0 = self.t0
        pkhdr.gv1 = self.gv1
        pkhdr.gv2 = self.gv2
        pkhdr.pv1 = self.pv1
        pkhdr.pv2 = self.pv2
        pkhdr.rmsd = self.rmsd_sem
        pkhdr.rmsr = self.rmsr_sem
        pkhdr.rmss = self.rmss_sem
        pkhdr.weight = self.weight_sem

        return pkhdr

    def clear_traces(self):
        """
        This clears out the 'data' part of the wavepacket, leaving
        only the header, and the various measures calculated from 
        the data.   Can use this after calling import_data and
        import_sem methods, which calculate the correlations, misfits 
        etc. if you don't need to retain the data, and the time 
        series calculated from it, to save space when manipulating
        the headers only 
        """

        if self._data_populated == True:
            del(self.data)
#            del(self.t)

        if self._sem_populated == True:
            del(self.sem_data)
#            del(self.ccorr)
#            del(self.ccorr_env)
#            del(self.ccorr_t)

        if self._data_env_populated:
            del(self.data_env)
#            del(self.data_iphase)

        if self._sem_env_populated:
            del(self.sem_env)
#            del(self.sem_iphase)

#        if self._data_env_populated and self._sem_env_populated:
#            del(self.env_ccorr)
#            del(self.env_ccorr_env)
#            del(self.env_ccorr_t)
            
        self._data_populated = False
        self._sem_populated = False
        self._data_env_populated = False
        self._sem_env_populated = False

    def calculate_ccorr_env_tmax(self):
        """ find time value (interpolating) of maximum of ccorr envelope """
        
        resamp_rate = 30 
        n_resamp = (len(self.ccorr_env) - 1)*resamp_rate + 1
        t_new = np.linspace(self.ccorr_t[0],self.ccorr_t[-1],n_resamp)
        ccenv_smooth = scipy.interpolate.spline(self.ccorr_t,self.ccorr_env,t_new)
        imax = np.argmax(ccenv_smooth)
        self.ccorr_env_tmax = t_new[imax]

    def __check_data(self):
        # check that real data are defined -- raise exception if not
        if not self._data_populated:
            raise Exception("Missing Data!")


    def __check_sem(self):  
        # Check that sem data are defined -- raise exception if not
        if not self._sem_populated:
            raise Exception("Missing SEM synthetic!")


    def __check_data_env(self):
        # check that data envelope is defined -- raise expection if not
        if not self._data_env_populated:
            raise Exception("Data envelope not calculated")


    def __check_sem_env(self):
        # check that sem envelope is defined -- raise exception if not
        if not self._sem_env_populated:
            raise Exception("SEM envelope not calculated")


    def __calculate_correlation_coeff(self):
        self.__check_data()
        self.__check_sem()

        # escalate warnings to errors -- calculate correlation coeff.
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.correlation = scipy.stats.pearsonr(self.data,self.sem_data)[0]


    def __calculate_env_correlation_coeff(self):
        self.__check_data_env()
        self.__check_sem_env()

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.env_correlation = scipy.stats.pearsonr(self.data_env,self.sem_env)[0]


    def __calculate_cross_correlation(self):
        self.__check_data()
        self.__check_sem()

        self.ccorr = np.correlate(self.data,self.sem_data,"full")
        self.ccorr_env = np.absolute(scipy.signal.hilbert(self.ccorr))
        self.ccorr_t = np.arange(-self.ndata+1,self.ndata)*self.smplintv


    def __calculate_data_env(self):
        hilb_data = scipy.signal.hilbert(self.data)
        self.data_env = np.absolute(hilb_data)
#        self.data_iphase = np.angle(hilb_data)
#        self.data_insfreq = pol.instantaneous_frequency(self.data, 1.)
#        nfft = 2**(np.ceil(np.log2(len(self.data)))+4)
#        nperseg = len(self.data)/3
#        noverlap = nperseg/2
        
#        [self.data_spec_f, self.data_spec_t, self.data_Sxx] = \
#                scipy.signal.spectrogram(self.data, 1./self.smplintv, \
#                nperseg=nperseg, noverlap=noverlap, nfft=nfft, mode='angle')
        
        self._data_env_populated = True


    def __calculate_sem_env(self):
        hilb_sem = scipy.signal.hilbert(self.sem_data)
        self.sem_env = np.absolute(hilb_sem)
#        self.sem_iphase = np.angle(hilb_sem)
#        self.sem_insfreq = pol.instantaneous_frequency(self.sem_data, 1.)
#        nfft = 2**(np.ceil(np.log2(len(self.sem_data)))+4)
#        nperseg = len(self.sem_data)/3
#        noverlap = nperseg/2
#        [self.sem_spec_f, self.sem_spec_t, self.sem_Sxx] = \
#                scipy.signal.spectrogram(self.sem_data, 1./self.smplintv, \
#                nperseg=nperseg, noverlap=noverlap, nfft=nfft, mode='angle')

#        self.sem_spect_psd=scipy.signal.spectrogram(self.sem_data, fs=1./self.smplintv,
#                                                    window=('tukey',.25), nperseg=256, noverlap=None,
#                                                    nfft=None, detrend='constant', return_onesided=True,
#                                                    scaling='density', axis=-1,mode='psd')
        self._sem_env_populated = True
  
    def __calculate_spec_correlation_coeff(self):
        self.__check_data()
        self.__check_sem()
        
        extval = 2
        
        dt = self.smplintv
    
        signal=self.data
        nsignal = 2**(int(np.log2(len(signal)))+extval)
        ft = np.fft.fft(signal,nsignal)
        
        signal=self.sem_data
        s_ft = np.fft.fft(signal,nsignal)
                    
        freqs = np.fft.fftfreq(nsignal, dt) * 1000.
        idx = np.argsort(freqs)
#        inds = np.where((freqs[idx] >= 1000./400.) & (freqs[idx] <= 1000./60.))[0] 
#        inds = np.where((freqs[idx] >= 4.) & (freqs[idx] <= 12.5))[0]
        inds = np.where((freqs[idx] >= 1000.*self.filt_w1/2./np.pi) & (freqs[idx] <= 1000*self.filt_w4/2./np.pi))[0]

        freqs = freqs[idx][inds]
        ft = ft[idx][inds]
        s_ft = s_ft[idx][inds]
        
        self.data_ft = ft
        self.sem_ft  = s_ft
        self.freqs = freqs
        
        self.sp_amp_correlation = scipy.stats.pearsonr(np.abs(ft),np.abs(s_ft))[0]
        self.sp_phase_correlation = scipy.stats.pearsonr(np.angle(ft),np.angle(s_ft))[0] 
        
        self.sp_amp_rmss = ((np.abs(ft))**2).sum()
        self.sp_amp_rmsd = ((np.abs(ft))**2).sum()
        
        self.sp_amp_rmss = np.sqrt(((np.abs(ft)-np.abs(s_ft))**2).sum()/self.sp_amp_rmss)
        self.sp_amp_rmsr = np.sqrt(((np.abs(ft)-np.abs(s_ft))**2).sum()/self.sp_amp_rmsd)
        self.sp_amp_rmsd = np.sqrt(self.sp_amp_rmsd/self.ndata)

        
        self.sp_phase_rmss = ((np.angle(ft))**2).sum()
        self.sp_phase_rmsd = ((np.angle(ft))**2).sum()
        
        self.sp_phase_rmss = np.sqrt(((np.angle(ft)-np.angle(s_ft))**2).sum()/self.sp_phase_rmss)
        self.sp_phase_rmsr = np.sqrt(((np.angle(ft)-np.angle(s_ft))**2).sum()/self.sp_phase_rmsd)
        self.sp_phase_rmsd = np.sqrt(self.sp_phase_rmsd/self.ndata)

        self.l2_sp_amp_misfit = ((np.abs(ft)-np.abs(s_ft))**2).sum()
        self.global_sp_amp_misfit = ((np.abs(ft)-np.abs(s_ft))**2).sum()*self.weight_data
        
#        self.l2_sp_amp_misfit = np.linalg.norm(np.abs(ft)-np.abs(s_ft))
        self.l2_sp_phase_misfit = np.linalg.norm(np.angle(ft)-np.angle(s_ft))
        
        self.data_smoothness_factor = np.std(np.diff(abs(ft)))/abs(np.mean(np.diff(abs(ft))))
        self.sem_smoothness_factor = np.std(np.diff(abs(s_ft)))/abs(np.mean(np.diff(abs(s_ft))))

#        # escalate warnings to errors -- calculate correlation coeff.
#        with warnings.catch_warnings():
#            warnings.simplefilter("error")
#            semspec=self.sem_Sxx.ravel()
#            dataspec=self.data_Sxx.ravel()
#            self.spec_correlation = np.sum(semspec*dataspec)/np.sqrt(np.sum(semspec*semspec)*np.sum(dataspec*dataspec))
            
    def __calculate_iphase_correlation_coeff(self):
        self.__check_data_env()
        self.__check_sem_env()

        # escalate warnings to errors -- calculate correlation coeff.
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.iphase_correlation = scipy.stats.pearsonr(self.data_iphase,self.sem_iphase)[0]

    def __calculate_insfreq_correlation_coeff(self):
        self.__check_data_env()
        self.__check_sem_env()

        # escalate warnings to errors -- calculate correlation coeff.
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.insfreq_correlation = scipy.stats.pearsonr(self.data_insfreq,self.sem_insfreq)[0]
            
    def __calculate_env_cross_correlation(self):
        self.__check_data_env()
        self.__check_sem_env()

        self.env_ccorr = np.correlate(self.data_env,self.sem_env,"full")
        self.env_ccorr_env = np.absolute(scipy.signal.hilbert(self.env_ccorr))
        self.env_ccorr_t = np.arange(-self.ndata+1,self.ndata)*self.smplintv


    def __calculate_l2_misfit(self):
        self.l2_misfit = np.linalg.norm(self.data-self.sem_data)
        self.data_sumsq = np.sqrt((self.data**2).sum())
        self.l2_misfit_normalized = self.l2_misfit / self.data_sumsq

    def __calculate_l1_misfit(self):
        self.l1_misfit = np.linalg.norm(self.data-self.sem_data,ord=1)
        self.data_sum = self.data.sum()
        self.l1_misfit_normalized = self.l1_misfit / self.data_sum

    def __calculate_env_l2_misfit(self):
        self.__check_data_env()
        self.__check_sem_env()
        
#        self.l2_env_misfit = np.linalg.norm(self.data_env-self.sem_env)
        self.l2_env_misfit = ((self.data_env-self.sem_env)**2).sum()
        self.data_env_sumsq = np.sqrt((self.data_env**2).sum())
        self.l2_env_misfit_normalized = self.l2_env_misfit / self.data_env_sumsq

    def __calculate_env_l1_misfit(self):
        self.l1_env_misfit = np.linalg.norm(self.data_env-self.sem_env,ord=1)
        self.data_env_sum = self.data_env.sum()
        self.l1_env_misfit_normalized = self.l1_env_misfit / self.data_env_sum
        
    def __calculate_gof(self):
        #(self.filt_w4/2./np.pi)
        self.pgf = misfit.pg(self.data,self.sem_data,self.smplintv, fmin=(self.filt_w1/2./np.pi), fmax=(self.filt_w4/2./np.pi), nf=100)
        self.egf = misfit.eg(self.data,self.sem_data,self.smplintv, fmin=(self.filt_w1/2./np.pi), fmax=(self.filt_w4/2./np.pi), nf=100)
        
#    def __calculate_global_misfit(self):
#        self.globalmisfit = (((self.data-self.sem_data)**2).sum())*self.weight_data
#        self.varred = np.sqrt(((self.data-self.sem_data)**2).sum())/np.sqrt(((self.data)**2).sum())
#        self.varredsem = np.sqrt(((self.data-self.sem_data)**2).sum())/np.sqrt(((self.sem_data)**2).sum())
#        self.sddiff = np.sum((self.data-self.sem_data)**2)
#        self.damp = np.sum((self.data)**2)

    def __calculate_global_misfit(self):
        self.globalmisfit = (((self.data-self.sem_data)**2).sum())*self.weight_data
        self.env_globalmisfit = (((self.data_env-self.sem_env)**2).sum())*self.weight_data
        
        self.u_rmss = ((self.sem_data)**2).sum()
        self.u_rmsd = ((self.data)**2).sum()
        
        self.u_rmss = np.sqrt(((self.data-self.sem_data)**2).sum()/self.u_rmss)
        self.u_rmsr = np.sqrt(((self.data-self.sem_data)**2).sum()/self.u_rmsd)
        self.u_rmsd = np.sqrt(self.u_rmsd/self.ndata)

        self.env_rmss = ((self.sem_env)**2).sum()
        self.env_rmsd = ((self.data_env)**2).sum()
        
        self.env_rmss = np.sqrt(((self.data_env-self.sem_env)**2).sum()/self.env_rmss)
        self.env_rmsr = np.sqrt(((self.data_env-self.sem_env)**2).sum()/self.env_rmsd)
        self.env_rmsd = np.sqrt(self.env_rmsd/self.ndata)

# ----------------------------------------------------------------------
# Wavepacket reading functions involving WavepacketData class
# ----------------------------------------------------------------------
def get_file_lists(wpDataDir,semDataDir=None):
    """given the directory containing wpd/wph files and that containing smd files,
    find all sets of matching wpd/wph/smd files and return lists containing
    their paths"""

    wpDirList = [ f for f in os.listdir(wpDataDir) 
                  if os.path.isfile(os.path.join(wpDataDir,f)) ]
    wphFiles = [ f for f in wpDirList if re.match(".*\.wph$",f) ]
    wphFiles.sort()

    wphMatches = []
    wpdMatches = []
    smdMatches = []

    for wphFile in wphFiles:
        wphPath = os.path.join(wpDataDir,wphFile)
        wpdFile = wphFile.replace(".wph",".wpd")
        wpdPath = os.path.join(wpDataDir,wpdFile)
        if semDataDir == None:
            if os.path.isfile(wpdPath):
                wphMatches.append(wphPath)
                wpdMatches.append(wpdPath)
        else:
            smdFile = wphFile.replace(".wph",".smd")
            smdPath = os.path.join(semDataDir,smdFile)
            if os.path.isfile(wpdPath) and os.path.isfile(smdPath):
                wphMatches.append(wphPath)
                wpdMatches.append(wpdPath)
                smdMatches.append(smdPath)

    if semDataDir == None:
        return wphMatches, wpdMatches, None
    else:
        return wphMatches, wpdMatches, smdMatches


def extract_wavepackets(wpkts_data,wpkts_sem,clear_traces=False):
    """
    combines data and sem wavepackets into a single list of 
    WavepacketData instances.   This function is called by 
    read_data on each event's wph/wpd/smd files
    
    If clear_traces == True then we remove the data from the 
    WavepacketData class, although we must first read it in 
    to calculate all the various correlations etc. that are done.
    """
    
    wpList = []

    # for each wavepacket, create WavepacketData object containing
    # all relevant header info, and both real and sem synthetic data
    for th_data,tw_data,th_sem,tw_sem in zip(wpkts_data.trace_hdrs,wpkts_data.traces,
                                             wpkts_sem.trace_hdrs,wpkts_sem.traces):
        for wpkt_data, wpkt_sem in zip(tw_data,tw_sem):
            new_wp = WavepacketData()
            new_wp.import_data(wpkts_data.event,th_data,wpkt_data)
            try:
                new_wp.import_sem(wpkt_sem)
                if clear_traces == True:
                    new_wp.clear_traces()
                wpList.append(new_wp)
            except RuntimeWarning:
                pass

    return wpList


def extract_wavepackets_nosem(wpkts_data,clear_traces=False):
    """
    Reads data into a single list of WavepacketData instances.   This function
    is called by read_data on each event's wph/wpd files
    
    If clear_traces == True then we remove the data from the WavepacketData
    class, although we must first read it in to calculate all the various
    correlations etc. that are done.
    """
    
    wpList = []

    # for each wavepacket, create WavepacketData object containing
    # all relevant header info, and both real and sem synthetic data
    for th_data,tw_data in zip(wpkts_data.trace_hdrs,wpkts_data.traces):
        for wpkt_data in tw_data:
            new_wp = WavepacketData()
            new_wp.import_data(wpkts_data.event,th_data,wpkt_data)
            if clear_traces == True:
                new_wp.clear_traces()
            wpList.append(new_wp)

    return wpList

def read_data(wphFiles,wpdFiles,smdFiles,clear_traces=False):
    """read all data into a list of WavepacketData instances"""

    wpList = []
    
    if smdFiles != None:
        for wphFile, wpdFile, smdFile in zip(wphFiles,wpdFiles,smdFiles):
            wpkts_data = read_UCB_data.wpData(wphFile,datafile=wpdFile)
            wpkts_data.read_data()
            wpkts_sem = read_UCB_data.wpData(wphFile,datafile=smdFile)
            wpkts_sem.read_data()
            wpList.extend(extract_wavepackets(wpkts_data,wpkts_sem,clear_traces))
    else:
        for wphFile, wpdFile in zip(wphFiles,wpdFiles):
            wpkts_data = read_UCB_data.wpData(wphFile,datafile=wpdFile)
            wpkts_data.read_data()
            wpList.extend(extract_wavepackets_nosem(wpkts_data,clear_traces))

    return wpList


def read_from_paths_pkl(wpDataDir=None,semDataDir=None,wp_pkl=None,clear_traces=False):
    
    # read pickle file if it exists and is specified
    if wp_pkl != None:
        if os.path.isfile(wp_pkl):
            wavepackets = cPickle.load(open(wp_pkl,'rb'))
            #print "Read in %d wavepackets from pickle file %s" % \
            #        (len(wavepackets),wp_pkl)
            return wavepackets
    
    # if pickle file not specified or doesn't exist, read directly
    if wpDataDir == None:
        print("Error: cannot find data to read!!")
        return None

    wphFiles, wpdFiles, smdFiles = get_file_lists(wpDataDir,semDataDir)
    wavepackets = read_data(wphFiles,wpdFiles,smdFiles,clear_traces)
    #print "Read in %d wavepackets" % len(wavepackets)
    if clear_traces:
        print("...retaining headers in memory but not traces")

    # save to pickle file if specified but it didn't exist
    if wp_pkl != None and (not os.path.exists(wp_pkl)):
        cPickle.dump(wavepackets,open(wp_pkl,'wb'),protocol=-1)
        
    return wavepackets


# ----------------------------------------------------------------------
# Functions working with lists of WavpacketData objects
# ----------------------------------------------------------------------
def get_values(wpList,var):
    """Loop over wpList and extract a list of the values
    of wp.<var>, in one to one correspondence with the 
    list wpList (c.f. get_unique_values)"""

    values = []
    for wp in wpList:
        try:
            val = getattr(wp,var)
        except AttributeError:
            print("Unknown class member: %s" % var)
            return None

        values.append(val)

    return values


def get_unique_values(wpList,var):
    """Loop over the list of WavepacketData instances in wpList and 
    extract the values of wp.<var>.  Return a list of the unique
    values and a list of their counts"""
    
    values = {}

    for wp in wpList:
        try:
            val = getattr(wp,var)
        except AttributeError:
            print("Unknown class member: %s" % var)
            return None, None

        try: 
            values[val] += 1
        except KeyError: 
            values[val] = 1

    return values.keys(), values.values()


def get_range(wpList,var):
    """ Get the minimum and maximum values of variable var in the 
    list of WavepacketData entries in wpList"""

    values = get_values(wpList,var)
    return min(values), max(values)

def clear_traces(wpList):
    """
    Clear time series from wavepackets, leaving just the headers.
    """

    for wp in wpList:
        wp.clear_traces()


def data_selector(wpList,test_type="EQ",**kwargs):
    """ each kwarg must have the keyword being a member of the class, and 
    the argument the value or list of values it can take.

    The result is a boolean 'and' between different kwargs (i.e. we can 
    select for having particular stations _and_ a particular phase), and 
    a boolean 'or' between list items for a given kwarg (i.e. can 
    choose to include wavepackets from several differenet stations
    
    Can make more complicated logical constructs by chaining 
    data_selector calls"""
        
    selectedList = []

    allowed_tests = ["EQ","LT","LE","GE","GT","NE"]
    if test_type not in allowed_tests:
        raise InputError("test_type does not take allowed value")
    
    # define equality testing according to function test_type argument
    if test_type == "EQ":
        def bool_test(val1,val2):
            return val1 == val2
    elif test_type == "LT":
        def bool_test(val1,val2):
            return val1 < val2
    elif test_type == "LE":
        def bool_test(val1,val2):
            return val1 <= val2
    elif test_type == "GE":
        def bool_test(val1,val2):
            return val1 >= val2
    elif test_type == "GT":
        def bool_test(val1,val2):
            return val1 > val2
    elif test_type == "NE":
        def bool_test(val1,val2):
            return val1 != val2
#    elif test_type == "IN":
#        def bool_test(val1,val2):
#            return val1 in val2

    for wp in wpList:
        # selection list keeps track of true/false for each
        # argument tested (kwargs).   Only accept the wavepacket
        # if test true for all kwargs ('and') - but 'or' between
        # alternative values in a list given as an argument.
        selections = []
        
        for key, item in kwargs.items():
            val = getattr(wp,key)
            selected = False

            if isinstance(item,basestring):
                if bool_test(val,item):
                    selected = True
            #elif hasattr(item,"__getitem__"):
            #    for it in item:
            #        if bool_test(val,it):
            #            selected = True 
            #            break
            #else:
            #    if bool_test(val,item):
            #        selected = True
            try:
                for it in item:
                    if bool_test(val,it):
                        selected = True
                        break
            except TypeError:
                if bool_test(val,item):
                    selected = True
    
            selections.append(selected)

        selections = np.array(selections)

        if selections.all():
            selectedList.append(wp)

    return selectedList


def corr_selector(wpList,mincorr,maxcorr):
    """ 
    select wavepackets with a correlation coefficient
    between data and SEM synthetics satisfying:
    mincorr < corr <= maxcorr
    """
    selected = data_selector(wpList,"GT",correlation=mincorr)
    selected = data_selector(selected,"LE",correlation=maxcorr)
    return selected


def path_selector(wpList,pathlist):
    """
    select given paths: pathlist is expected to be a numpy structured
    array with fields 'event', 'stn' and 'phase' present
    """

    selectedList = []

    for path in pathlist:
        for wp in wpList:
            if ((path['event'] == wp.event_id) and
                (path['stn'] == wp.stn) and
                (path['phase'] == wp.phase) and
                (wp.weight_data > 0)):
                selectedList.append(wp)

    return selectedList


def sort_wavepackets(wpList,var,reverse=False):
    """
    return wpList sorted based on variable var (must be one
    of the attributes of the WavepacketData class)
    reverse == True implements reverse sorting
    """
    
    wpListSorted = sorted(wpList,
                          key=lambda x: getattr(x,var),
                          reverse=reverse)

    return wpListSorted

def sort_remove_outliers_wpkts(wpList,var,pc_bot,pc_top):
    """
    sort wpList according to variable var (using sort_wavepackets
    function), and then remove the top pc_top % and bottom pc_bot %
    of outliers.
    """

    wpListSorted = sort_wavepackets(wpList,var)
    
    n_list = len(wpListSorted)
    n_top = int(math.floor(n_list*pc_top/100.))
    n_bot = int(math.floor(n_list*pc_bot/100.))

    print("Variable: %s, outliers removed from top: %d, bottom: %d" % \
            (var,n_top,n_bot))

    return wpListSorted[n_bot:(n_list-n_top)]


def write_paths(wpList,outfile):
    """
    write path details (event, stn, event_lat, event_lon, stn_lat, stn_lon)
    to file, for all wavepackets in wpList
    """
    
    outf = open(outfile,"w")
    for wp in wpList:
        outf.write("%s %s %f %f %f %f\n" % (wp.event_id,wp.stn,
                                            wp.event_lat,wp.event_lon,
                                            wp.stn_lat,wp.stn_lon))
    outf.close()
    
def write_paths_phases(wpList,outfile):
    """
    write path details (event, stn, event_lat, event_lon, stn_lat, stn_lon)
    to file, for all wavepackets in wpList
    """
    
    outf = open(outfile,"w")
    for wp in wpList:
        outf.write("%s %s %d\n" % (wp.event_id,wp.stn,wp.phase))
    outf.close()

def write_paths_with_corr(wpList,outfile):
    """
    write path details (event, stn, event_lat, event_lon, stn_lat, stn_lon)
    to file, for all wavepackets in wpList
    """
    
    outf = open(outfile,"w")
    outf.write("event \t stn \t\t evlat \t\t evlon \t\t stlat \t\t stlon \t\t corr \t\t env_corr \t pgf \t egf \t data-weight \t sem-weight\n")
            
    for wp in wpList:
            outf.write("%s %s \t\t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n" % (wp.event_id,wp.stn,
                                        wp.event_lat,wp.event_lon,
                                        wp.stn_lat,wp.stn_lon,
                                        wp.correlation,wp.env_correlation,
                                        wp.weight_data, wp.weight_sem))
#        outf.write("%s %s \t\t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n" % (wp.event_id,wp.stn,
#                                            wp.event_lat,wp.event_lon,
#                                            wp.stn_lat,wp.stn_lon,
#                                            wp.correlation,wp.env_correlation,wp.pgf,
#                                            wp.egf,
#                                            wp.weight_data, wp.weight_sem))
    outf.close()

def write_paths_with_corr_phase(wpList,outfile):
    """
    write path details (event, stn, event_lat, event_lon, stn_lat, stn_lon)
    to file, for all wavepackets in wpList
    """
    
    outf = open(outfile,"w")
    outf.write("event \t\t stn \t netwk \t chnnl \t phase \t evlat \t\t evlon \t\t stlat \t\t stlon \t\t corr \t\t env_corr \t phase_corr \t insfreq_corr \t data-weight \t sem-weight\n")
            
    for wp in wpList:
        outf.write("%s \t %s \t %s \t %s \t %d \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n" % (wp.event_id,wp.stn,wp.netwk,wp.chnnl,wp.phase,
                                            wp.event_lat,wp.event_lon,
                                            wp.stn_lat,wp.stn_lon,
                                            wp.correlation,wp.env_correlation,wp.iphase_correlation,
                                            wp.insfreq_correlation,
                                            wp.weight_data, wp.weight_sem))
    outf.close()
 
def write_paths_with_gof_misfit(wpList,outfile):
    """
    write path details (event, stn, event_lat, event_lon, stn_lat, stn_lon)
    to file, for all wavepackets in wpList
    """
    
    outf = open(outfile,"w")
    outf.write("event \t stn \t comp \t chnnl \t phase \t evlat \t evlon \t stlat \t stlon \t corr \t env_corr \t pgf \t egf \t globalmisfit \t l2_misfit \t l2_misfit_normalized \t l2_env_misfit \t l2_env_misfit_normalized \t rmsr_data \t rmsd_data \t rmss_data \t rmsr_sem \t rmsd_sem \t rmss_sem \t weight_data \t weight_sem \t l1_misfit \t l1_misfit_normalized \t l1_env_misfit \t l1_env_misfit_normalized\n")
            
    for wp in wpList:
        outf.write("%s \t %s \t %c \t %s \t %d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e\n" 
                                         % (wp.event_id,wp.stn,wp.comp,wp.chnnl,wp.phase,
                                            wp.event_lat,wp.event_lon,
                                            wp.stn_lat,wp.stn_lon,
                                            wp.correlation,wp.env_correlation,wp.pgf,
                                            wp.egf, wp.globalmisfit,
                                            wp.l2_misfit, wp.l2_misfit_normalized,
                                            wp.l2_env_misfit, wp.l2_env_misfit_normalized,
                                            wp.rmsr_data, wp.rmsd_data, wp.rmss_data,
                                            wp.rmsr_sem, wp.rmsd_sem, wp.rmss_sem,
                                            wp.weight_data, wp.weight_sem,
                                            wp.l1_misfit, wp.l1_misfit_normalized,
                                            wp.l1_env_misfit, wp.l1_env_misfit_normalized))
    outf.close()    
    
def write_paths_with_rmsr_rmss(wpList,outfile):
    """
    write path details (event, stn, event_lat, event_lon, stn_lat, stn_lon)
    to file, for all wavepackets in wpList
    """
    
    outf = open(outfile,"w")
    outf.write("event \t stn \t\t evlat \t\t evlon \t\t stlat \t\t stlon\
\t\t phase \t\t rmsr \t\t rmss \t\t env_rmsr \t\t env_rmss \t\t sp_amp_rmsr \t\t sp_amp_rmss \t\t sp_phase_rmsr \t\t sp_phase_rmss\
\t\t correlation \t\t env_correlation \t\t sp_amp_correlation \t\t sp_phase_correlation\
\t\t weight\n")
            
    for wp in wpList:
            outf.write("%s %s \t\t %f \t %f \t %f \t %f \t %i \t %f \t %f \t %f \t %f \t \
%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n" % (wp.event_id,wp.stn,
                                        wp.event_lat,wp.event_lon,
                                        wp.stn_lat,wp.stn_lon,wp.phase,
                                        wp.u_rmsr,wp.u_rmss,
                                        wp.env_rmsr, wp.env_rmss,
                                        wp.sp_amp_rmsr, wp.sp_amp_rmss,
                                        wp.sp_phase_rmsr, wp.sp_phase_rmss,
                                        wp.correlation, wp.env_correlation,
                                        wp.sp_amp_correlation, wp.sp_phase_correlation,
                                        wp.weight_data))
#        outf.write("%s %s \t\t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n" % (wp.event_id,wp.stn,
#                                            wp.event_lat,wp.event_lon,
#                                            wp.stn_lat,wp.stn_lon,
#                                            wp.correlation,wp.env_correlation,wp.pgf,
#                                            wp.egf,
#                                            wp.weight_data, wp.weight_sem))
    outf.close()
    
def run_calculate_ccorr_env_tmax(wplist):

    i = 0
    nwp = len(wplist)

    for wp in wplist:
        wp.calculate_ccorr_env_tmax()
        i += 1
        if (i%100) == 0:
            print("Calculating time of cross-corr env max: %d of %d" % (i,nwp))
            
#
#def return_corr_list(wpList):
#    """
#    calculate the correlation between the SEM and acal synthetics
#    """
#    
#    for wp in wpList:
#        if wp._sem_populated == False:
#            print "Error: no SEM synthetics present"
#            return
#
#    # calculate correlation and store results along with path details
#    # in list
#    corr_list = []
#    
##    for wp in wpList:
##        corr_list.append( (wp.event_id,wp.stn,wp.phase,wp.correlation,wp.iphase_correlation,wp.env_correlation,wp.insfreq_correlation,wp.spec_correlation) )
##
##    # convert the list to a structured numpy array
##    corr_arr_dtype = np.dtype([('event','a8'),
##                               ('stn','a4'),
##                               ('phase',np.int32),
##                               ('correlation',np.float32),
##                               ('iphase_correlation',np.float32),
##                               ('env_correlation',np.float32),
##                               ('insfreq_correlation',np.float32),
##                               ('spec_correlation',np.float32)])
#    
#    for wp in wpList:
#        corr_list.append( (wp.event_id,wp.stn,wp.phase,wp.correlation,wp.env_correlation,wp.pgf,wp.egf) )
#
#    # convert the list to a structured numpy array
#    corr_arr_dtype = np.dtype([('event','a8'),
#                               ('stn','a4'),
#                               ('phase',np.int32),
#                               ('correlation',np.float32),
#                               ('env_correlation',np.float32),
#                               ('pgf',np.float32),
#                               ('egf',np.float32)])
#
#    corr_arr = np.zeros(len(corr_list),dtype=corr_arr_dtype)
#    for i,elem in enumerate(corr_list):
#        corr_arr[i] = elem
#
#    # sort by correlation
##    corr_arr.sort(order='corr')
##    corr_arr = corr_arr[::-1]
#    return corr_arr    
    
    
def return_corr_list(wpList):
    """
    calculate the correlation between the SEM and acal synthetics
    """

    for wp in wpList:
        if wp._sem_populated == False:
            print("Error: no SEM synthetics present")
            return

    # calculate correlation and store results along with path details
    # in list
    corr_list = []

#    for wp in wpList:
#        corr_list.append( (wp.event_id,wp.stn,wp.phase,wp.correlation,wp.iphase_correlation,wp.env_correlation,wp.insfreq_correlation,wp.spec_correlation) )
#
#    # convert the list to a structured numpy array
#    corr_arr_dtype = np.dtype([('event','a8'),
#                               ('stn','a4'),
#                               ('phase',np.int32),
#                               ('correlation',np.float32),
#                               ('iphase_correlation',np.float32),
#                               ('env_correlation',np.float32),
#                               ('insfreq_correlation',np.float32),
#                               ('spec_correlation',np.float32)])

    for wp in wpList:
        corr_list.append( (wp.event_id,wp.stn,wp.phase,wp.correlation,wp.env_correlation,wp.pgf,wp.egf,wp.globalmisfit,wp.varred,wp.varredsem,wp.sddiff,wp.damp) )

    # convert the list to a structured numpy array
    corr_arr_dtype = np.dtype([('event','a8'),
                               ('stn','a4'),
                               ('phase',np.int32),
                               ('correlation',np.float32),
                               ('env_correlation',np.float32),
                               ('pgf',np.float32),
                               ('egf',np.float32),
                               ('gmisfit',np.float32),
                               ('varred',np.float32),
                               ('varredsem',np.float32),
                               ('sddiff',np.float32),
                               ('damp',np.float32)])

    corr_arr = np.zeros(len(corr_list),dtype=corr_arr_dtype)
    for i,elem in enumerate(corr_list):
        corr_arr[i] = elem

    # sort by correlation
#    corr_arr.sort(order='corr')
#    corr_arr = corr_arr[::-1]
    return corr_arr

