'''
    Scan for drugs
'''
from ReadML import HighestPeaks, SelectPeaks
import pymzml
from pprint import pprint
import Tkinter, tkFileDialog

def GetRTByMZList(filename, mz_list, tolerance = 0.2):
    # select the retention time with most peaks
    # ignore retention time with less than 50% peaks of mz_list length
    run = pymzml.run.Reader(filename, noiseThreshold = 100)
    rt_peaks = dict() # peaks at special rt time
    max_len  = 0
    for spec in run:
        try:
            rt_time = spec["scan time"]
        except:
            continue
        rt_peak_list = []
        for eachmz in mz_list:
            eachrange = (eachmz - tolerance, eachmz + tolerance)
            try:
                mz, intensity = HighestPeaks(SelectPeaks(spec.peaks, eachrange))
            except Exception as e:
                #print e.message
                continue
            if intensity > 0:
                rt_peak_list.append([mz, intensity])
        if len(rt_peak_list) > 0.5 * len(mz_list):
            if len(rt_peak_list) > max_len:
                max_len = len(rt_peak_list)
            rt_peaks[rt_time] = rt_peak_list
    rt_peaks_reduce = dict()
    for mz in rt_peaks:
        if len(rt_peaks[mz]) == max_len:
            rt_peaks_reduce[mz] = rt_peaks[mz]
    return rt_peaks_reduce

def main():
    root = Tkinter.Tk()
    root.withdraw()
    inputfile = tkFileDialog.askopenfilename()
    peak_reduce = GetRTByMZList(inputfile, [423.3,405.2,296.2,268.2,235.2,206.1,171.1,143.1,115.1])
    pprint(peak_reduce)
    print sorted(peak_reduce.keys())

if __name__ == "__main__":
    main()
