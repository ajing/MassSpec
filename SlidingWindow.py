'''
    Sliding window algorithm for tRNA unique mass
    function need attention: removeNoise, deconvolute_peaks
'''

from ExtractChrom import ExtractSpec
from Common import InRange
from GetPeaks import SelectPeaks, HighestPeaks
from RNAMassGen import GetUniqueMass
import pymzml
from pprint import pprint

def AllMassPossibility(a_mass):
    return [a_mass, a_mass / 2, a_mass / 3, a_mass / 4]

def GetMaxPeakInSpecs(mz_list, specs, run, tolerance = 0.2):
    # get the maximun peaks in a collection of specs
    max_int_dict = dict()
    for eachmz in mz_list:
        max_int_dict[eachmz] = {"max_int": 0, "max_mz": None, "max_time": None, "max_id": None}
    max_intensity = 0
    for spec_basic in specs:
        idx  = spec_basic.index
        spec = run[idx]
        for eachmz in mz_list:
            eachrange = (eachmz - tolerance, eachmz + tolerance)
            try:
                mz, intensity = HighestPeaks(SelectPeaks(spec.peaks, eachrange))
            except Exception as e:
                #print e.message
                continue
            #print spec_basic
            #print eachrange
            #print "scan time:", spec["scan time"]
            max_intensity = max_int_dict[eachmz]["max_int"]
            if intensity > max_intensity:
                #print "intensity:", intensity
                max_int_dict[eachmz] = {"max_int": intensity, "max_mz": mz, "max_time": spec["scan time"], "max_id": spec["id"]}
    return max_int_dict

class MostPeaks:
    # Rule to pick the best match in mass spec with a list of mass
    # Basically, find the one with most peaks. The way to break tie is sum of intensity
    def __init__(self):
        self._max_num_peaks = 0
        self._max_intensity = 0
        self._cur_max_dict  = None

    def compare(self, new_int_dict):
        new_num_peaks = len([x for x in new_int_dict if new_int_dict[x]["max_int"] > 0])
        if new_num_peaks == 0:
            return
        new_intensity = sum([new_int_dict[x]["max_int"] for x in new_int_dict.keys()])
        if new_num_peaks > self._max_num_peaks or (new_num_peaks == self._max_num_peaks and new_intensity > self._max_intensity):
            #print new_int_dict
            self._max_num_peaks = new_num_peaks
            self._max_intensity = new_intensity
            self._cur_max_dict  = new_int_dict


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def SlidingWindow(masslist, run, exspec, rtrange = None, s_win = 5):
    # s_win: sliding window size, which is in minute
    rule   = MostPeaks()
    for rt_time in drange(exspec.start_time, exspec.end_time, 0.1):
        # ignore the spec out of rtrange
        if (not rtrange is None) and (rt_time < min(rtrange) or rt_time > max(rtrange)):
            continue
        print rt_time
        specs  = exspec.extractWithTimeRange(rt_time, rt_time + s_win)
        max_int_dict = GetMaxPeakInSpecs(masslist, specs, run)
        rule.compare(max_int_dict)
    print rule._cur_max_dict

def main():
    ms_file   = "./E165ug.mzML"
    # data structure for mass_dict
    # {'tRNA | Ile | CAU | Escherichia coli | prokaryotic cytosol': [['CCCUU4AGp', 2541.3703], ['ACU}AU6APCGp', 3786.716000000001]]}
    mass_dict = GetUniqueMass("./Data/tRNAseq.txt")
    exspec = ExtractSpec(ms_file)
    mass_all_possible = dict()
    for rna in mass_dict:
        run = pymzml.run.Reader(ms_file, noiseThreshold = 100)
        mass_list = []
        mass_all_possible[rna] = dict()
        for seg_mass in mass_dict[rna]:
            all_mass = AllMassPossibility(seg_mass[1])
            mass_all_possible[rna][seg_mass[0]] = all_mass
            mass_list += all_mass
        #print mass_dict[rna]
        #print mass_list
        print rna
        SlidingWindow(mass_list, run, exspec)


if __name__ == "__main__":
    main()
