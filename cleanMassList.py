'''
    Clean mass list file
'''

import csv

def CleanMassList(infile):
    outfile = infile + "_mod"
    outobj  = csv.writer(open(outfile, "w"), quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for line in open(infile):
        content = line.strip().split()
        if len(content) == 6:
            content.insert(2, "")
        outobj.writerow([content[0], content[3], content[6]])

if __name__ == "__main__":
    CleanMassList("modified_bases_list.txt")
