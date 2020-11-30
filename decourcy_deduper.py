Use this bash command to sort the file by chromosome number, prior to running the deduper python script
samtools sort -o test_sorted.sam test.sam

======================================

# Deduper script
#!/usr/bin/env python

import re
import argparse

parser = argparse.ArgumentParser(description="'-h' takes the following args:")
parser.add_argument("-f", "--in_file", help="Enter input file - absolute file path and file name", required=True)
parser.add_argument("-o", "--out_file_prefix", help="Enter output file prefix", required=True)
parser.add_argument("-p", "--paired", help="Is file paired end (not single end): T/F", required=False)
parser.add_argument("-u", "--umi", help="Enter UMI list file - absolute file path and name", required=False)

args = parser.parse_args()

r_in = args.in_file
r_out_pre = args.out_file_prefix
r_paired = args.paired
r_umi = args.umi


def cigar_adj(cigar):
    chunked = (re.findall("[0-9]+[A-Z]{1}", cigar))

    if chunked[0][-1] == "M":
        return(0)

    # Go through and subtract position numbers if they are from S
    elif chunked[0][-1] == "S":
        return(int(chunked[0][:-1]))


# Takes an integer - the bit flag - and returns T/F based on whether
# the "SEQ being reverse complemented" portion of the flag is true
# in the boolean variable isrev (IS REVerse complemented)
def bit_flag_interpreter(an_int):
    isrev = None
    if ((an_int & 16) == 16):
        isrev = True
    else:
        isrev = False
    
    return(isrev)


# setting empty variables for use later in the program
add_list = []
prev_chrom = ''
check_set = set()



with open(r_in, 'rt') as f, open("%s_removed_reads_deduped.sam" % r_out_pre, 'w') as removed, open("%s_kept_reads_deduped.sam" % r_out_pre, 'w') as kept:
    for line in f:
        if "@" in line.split()[0]:
            kept.write(line)
            removed.write(line)
        if "@" not in line.split()[0]:
            chromosome = (line.split()[2])

            if chromosome != prev_chrom:
                check_set = set()
            # if the previously used chromosome isn't the same as this one, this clears the 
            # set being used to keep track of PCR duplicates
            # if they are different chromosome numbers, this clears the list.

            umi = (line.split()[0].split(":")[-1])
            bit_flag = int(line.split()[1])            
            raw_location = int(line.split()[3])
            cigar = (line.split()[5])


            # Converts the raw location to the true location based on the cigar string
            # and possible soft-clippings if needed - or passes the raw -> true location if not
            true_location = raw_location + cigar_adj(cigar)

            # revc is boolean variable output T/F:
            # True = read is reverse complement
            revc = bit_flag_interpreter(bit_flag)

            # wlist - W(orking) list - holds the info needed to check reads for being duplicates
            wlist = [chromosome, umi, true_location, revc]
            

            if tuple(wlist[1:4]) in check_set:
                # write to removed duplicates file            
                removed.write(line)

            elif tuple(wlist[1:4]) not in check_set:
                # write to "keep" file
                # add to check_set
                check_set.add(tuple(wlist[1:4]))
                kept.write(line)

            prev_chrom = chromosome
            # ^ Sets the previous chromsome so that the program can check if a new
            # chromosome's PCR reads are being analyzed

            # close files
f.close()
removed.close()
kept.close()

