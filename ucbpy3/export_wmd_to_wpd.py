#!/usr/bin/env python

#path2ucbpy = '/home/gcl/BR/fmunch/ucbpy3' 
## Load ucbpy3 (Python 3 version of ucbpy) 
#import sys 
#sys.path.insert(0, path2ucbpy) 

import re
import argparse
from os.path import splitext
from WPDataRW import WPData

encoding = 'ISO-8859-1'#utf-8
#encoding = 'utf-8'

packhdr_exts = {
    'smd': 'packhdrbr_st', 
    'wpd': 'packhdrbr_st', 
    'wmd': 'packhdrmp_st'}
tracehdr_exts = {
    'wph': 'tracehdr_st',
    'wmh': 'tracehdr_st',
    'wpH': 'tracehdrH_st'}

def export_to_wpd(events, all_packets, path):
    for event in events:
        fh_in  = '%s/%s.wmh' % (path, event)
        fd_in  = '%s/%s.wmd' % (path, event)
        fh_out = '%s/%s.wph' % (path, event)
        fd_out = '%s/%s.wpd' % (path, event)
        wpd_in = WPData(fh_in, fd_in,
            tracehdr_type=tracehdr_exts['wmh'],
            packhdr_type=packhdr_exts['wmd'])
        wpd_out = WPData(fh_out, fd_out,
            tracehdr_type=tracehdr_exts['wph'],
            packhdr_type=packhdr_exts['wpd'])
        if wpd_in.open('r') and wpd_out.open('w'):
            title, traces = wpd_in.read_headers()
            traces_out = []
            ntr_kept = 0
            nwp_kept = 0
            for trace in traces:
                packets_out = []
                for packhdr, data in wpd_in.read_data(trace):
                    if all_packets or int(packhdr['keep']) == 1:
                        packets_out.append((packhdr,data))
                        nwp_kept += 1
                if packets_out:
                    traces_out.append((trace, packets_out))
                    ntr_kept += 1
            if traces_out:
                wpd_out.write_headers_and_data(title, traces_out)
                print('Kept %i wavepackets (%i traces) for event %s' % (
                    nwp_kept, ntr_kept, event) )
            else:
                print('Waring: No traces kept for event %s' % (event) )
            wpd_out.close()
            wpd_in.close()
        else:
            print('Warning: could not open files for event %s - skipping ...' % (event))

def main():
    # set up parser for command-line args
    parser = argparse.ArgumentParser(description='Export wmh/wmd - type wavepackets as wph/wpd')
    parser.add_argument('-a', '--all-packets', action='store_true',
        help='export all wavepackets (including those _not_ marked "keep")')
    parser.add_argument('-p', '--path', type=str, default='.',
        help='path for reading / writing wavepackets')
    parser.add_argument('events', metavar='event', type=str, nargs='+',
        help='event codes to process')
    arg = parser.parse_args()
    export_to_wpd(arg.events, arg.all_packets, arg.path)

if __name__ == '__main__':
    main()
