'''
    Generate tRNA mass and corresponding sequence
'''
import re
import csv
from itertools import groupby

MASS_FILE = "./Data/modified_bases_list.txt_mod"



def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def MassListParser(infile):
    # take mass list and return dictionary with seq as key and mass as value
    mass_dict = dict()
    input_file = csv.DictReader(open(infile), quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for row in input_file:
        for eachkey in row:
            if is_number(row[eachkey]) and eachkey != "symbol":
                row[eachkey] = float(row[eachkey])
        mass_dict[row["symbol"]] = row
    return mass_dict

def GetMassForBase(base, mass_dict):
    if len(base) == 1:
        return mass_dict[base]["N"]+1.0079
    elif base.endswith("p"):
        return mass_dict[base[0]]["N"]+80.998

def GetMassForSeq(seq, mass_dict):
    total = 0
    seg = re.findall(".[p]?", seq)
    for each in seg:
        #print each,GetMassForBase(each, mass_dict)
        total += GetMassForBase(each, mass_dict)
    total += (len(seq.replace("p","")) - 1) * 61.965
    return total

def FivePrimePho(seq):
    # five prime phosphate
    return seq + "p"

def ThreePrimePho(seq):
    # three prime phosphate
    return "p" + seq

def RNAGen(seq):
    # take tRNA sequence as input and generate all unique segment aftersplicing
    # RNase T1 digest
    seg = seq.split("G")
    #seg = re.findall(seq, "[^G]*G")
    # the last element in list is not spliced by digest protein
    seg_mod = [x + "Gp" for x in seg[:-1]]
    seg_mod.append(seg[-1])
    final_seg = list(set(seg_mod))
    final_seg.sort(key=len)
    return final_seg

def MassList(seq, mass_dict):
    seq_list = RNAGen(seq)
    masslist = []
    for each in seq_list:
        masslist.append([each, GetMassForSeq(each, mass_dict)])
    return masslist

# FastA file parser
def FastAIter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        seq = seq.replace("-", "")
        yield header, seq

def UniqueMass(mass_list, tolerance):
    mass_sort = sorted(mass_list)
    uniq_elem = []
    for i in range(len(mass_sort) - 1):
        if i >= 1:
            if not (mass_sort[i + 1] - mass_sort[i] < tolerance or mass_sort[i] - mass_sort[i - 1] < tolerance):
                return value
    return False

def UniqueSegMass(seg_mass_list, tolerance = 0.5):
    seg_mass = dict()
    mass_list = [ mass for segname, mass in mass_list]
    for seg_name, mass in mass_list:
        for key_mass in seg_mass:
        if not mass in seg_mass:
            seg_mass[mass] = seg_name
        else:
            seg_mass[mass] = None
    seg_mass_mod = dict()
    for key, value in seg_mass.iteritems():
        if not value is None:
            seg_mass_mod[key] = value
    return seg_mass_mod.keys()

def Test():
    massdict =  MassListParser(MASS_FILE)
    #GetMassForSeq("GGGp", massdict)
    #print GetMassForSeq("CD", massdict)
    #assert(RNAGen("GGGGCUAUAGCUCAGCD") == ["CD", "Gp" , "CUCAGp", "CUAUAGp"])
    #MassList("GGGGCUAUAGCUCAGCD", massdict)
    each_trna = dict()
    for header, seq in FastAIter("./Data/tRNAseq.txt"):
        #print header, "seq", seq
        each_trna[header] = dict()
        each_trna[header]["seq"] = seq
        each_trna[header]["mass_list"] = MassList(seq, massdict)
    seg_mass = []
    for key, value in each_trna.iteritems():
        seg_mass += value["mass_list"]
    print seg_mass, "yes"
    print "yes"

if __name__ == "__main__":
    Test()
    #RNAGen("GGGGCUAUAGCUCAGCDGGGAGAGCGCCUGCUUVGCACGCAGGAG")
