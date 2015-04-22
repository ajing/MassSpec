'''
    Generate tRNA mass and corresponding sequence
'''
import re
import csv
import Tkinter, tkFileDialog
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

def UniqueMass(mass_list, tolerance = 0.5):
    mass_sort = sorted(mass_list)
    uniq_elem = []
    i = 0
    while i < len(mass_sort):
        if i < len(mass_sort) - 1 and mass_sort[i + 1] - mass_sort[i] <= tolerance:
            i += 2
        elif i > 0 and mass_sort[i] - mass_sort[i - 1] <= tolerance:
            i += 1
        else:
            uniq_elem.append(mass_sort[i])
            i += 1
    #print "final elements:", uniq_elem
    return uniq_elem

def UniqueRNAMass(rna_dict, uniq_m_list):
    rna_uniq_mass = dict()
    for header in rna_dict:
        rna_uniq_mass[header] = []
        for seg, mass in rna_dict[header]["mass_list"]:
            if mass in uniq_m_list:
                rna_uniq_mass[header].append([seg, mass])
    return rna_uniq_mass

def PrintNiceRNAMass(rna_dict):
    for header in rna_dict:
        print header
        print "Unique Mass: " + "\t".join([",".join(map(str, x)) for x in rna_dict[header]])

def Test():
    #GetMassForSeq("GGGp", massdict)
    #print GetMassForSeq("CD", massdict)
    #assert(RNAGen("GGGGCUAUAGCUCAGCD") == ["CD", "Gp" , "CUCAGp", "CUAUAGp"])
    #MassList("GGGGCUAUAGCUCAGCD", massdict)
    assert UniqueMass([1,5,6,9], 1) == [1,9], UniqueMass([1,5,6,9], 1)
    assert UniqueMass([1,2,6,9], 1) == [6,9], UniqueMass([1,2,6,9], 1)
    assert UniqueMass([1,2,3,6,9], 1) == [6,9], UniqueMass([1,2,3,6,9], 1)
    assert UniqueMass([1,2,4,8,9], 1) == [4], UniqueMass([1,2,4,8,9], 1)

def GetUniqueMass(inputfile):
    massdict =  MassListParser(MASS_FILE)
    each_trna = dict()
    for header, seq in FastAIter(inputfile):
        #print header, "seq", seq
        each_trna[header] = dict()
        each_trna[header]["seq"] = seq
        each_trna[header]["mass_list"] = MassList(seq, massdict)
    seg_list = []
    for key, value in each_trna.iteritems():
        seg_list += value["mass_list"]
    mass_list = [x[1] for x in seg_list]
    uniq_m_list    = UniqueMass(mass_list)
    uniq_rna_mass  = UniqueRNAMass(each_trna, uniq_m_list)
    return uniq_rna_mass

def TestSeg(inputfile):
    uniq_rna_mass = GetUniqueMass(inputfile)
    PrintNiceRNAMass(uniq_rna_mass)

if __name__ == "__main__":
    root = Tkinter.Tk()
    root.withdraw()
    inputfile = tkFileDialog.askopenfilename()
    TestSeg(inputfile)
    #RNAGen("GGGGCUAUAGCUCAGCDGGGAGAGCGCCUGCUUVGCACGCAGGAG")
