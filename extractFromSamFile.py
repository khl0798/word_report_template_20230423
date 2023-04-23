#!coding: utf-8
import pysam
import os,re
# from pysam import utils

sam_path="/mount/labdatastorage001/input/konghl/data/plant/MiS027/results_os_new/sam_file"
samfile2="MiS027_ANS01WK_27x03x2023_0_0031xR2_0_LTA_0_ANS01WKx0031_2_0_MiS027_F_M235_500_BE1017_CTGTGAGATG_GTGTTATTAAGTTGTCTAAGCGTC_OryzasativaL_TDNA_CCTAACTGCTGTGCCACT_CTGCTACGATGA_TAGTCATCGTAGCAGCCTAAC_.fastq.trimmed2.sam"
samfile1="MiS027_ANS01WK_27x03x2023_0_0031xR1_0_LTA_0_ANS01WKx0031_1_0_MiS027_F_M235_500_BE1016_CTGTGAGATG_GTGTTATTAAGTTGTCTAAGCGTC_OryzasativaL_TDNA_CCTAACTGCTGTGCCACT_AGTCACTGAGAG_TAGCTCTCAGTGACTCCTAAC_.fastq.trimmed2.sam"
samfile3="MiS027_ANS01WK_27x03x2023_0_0031xR3_0_LTA_0_ANS01WKx0031_3_0_MiS027_F_M235_500_BE1018_CTGTGAGATG_GTGTTATTAAGTTGTCTAAGCGTC_OryzasativaL_TDNA_CCTAACTGCTGTGCCACT_ACGCTGAGAGAT_TAGATCTCTCAGCGTCCTAAC_.fastq.trimmed2.sam"
samfile_list=[samfile1,samfile2,samfile3]


class CustomAlignedSegment:
    def __init__(self):
        self.query_name = None
        self.flag = None
        self.reference_name = None
        self.reference_start = None
        self.mapping_quality = None
        self.cigarstring = None
        self.next_reference_name = None
        self.next_reference_start = None
        self.template_length = None
        self.query_sequence = None
        self.query_qualities = None

def parse_sam_line(sam_line):
    fields = sam_line.strip().split("\t")
    if len(fields) < 11:
        raise ValueError("Invalid SAM line")
    read = CustomAlignedSegment()
    read.query_name = fields[0]
    read.flag = int(fields[1])
    read.reference_name = fields[2]
    read.reference_start = int(fields[3]) - 1
    read.mapping_quality = int(fields[4])
    read.cigarstring = fields[5]
    read.next_reference_name = fields[6]
    read.next_reference_start = int(fields[7]) - 1
    read.template_length = int(fields[8])
    read.query_sequence = fields[9]
    # read.query_qualities = pysam.qualitystring_to_array(fields[10])
    read.query_qualities=str(fields[10])
    return read

def parse_cigar(cigarstring):
    cigar_operations = ('M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X')
    cigar_re = re.compile(r"(\d+)([{}])".format("".join(cigar_operations)))
    return [(int(length), op) for length, op in cigar_re.findall(cigarstring)]


def calculate_reference_end(read):
    if read.cigarstring is None:
        return None
    cigar_tuples = parse_cigar(read.cigarstring)
    reference_end = read.reference_start
    for length, op in cigar_tuples:
        if op in ('M', 'D', 'N', '=', 'X'):  # M, D, N, =, X
            reference_end += int(length)
    return reference_end - 1


# filename = os.path.join(sam_path,samfile)

# outfile="test1.sam"
tmp1=[]
for s in samfile_list:
    filename=os.path.join(sam_path, s)
    print(s)
    with open(filename, "r") as samfile:
        for line in samfile:
            if line.startswith("@"):
                continue  # Ignore header lines
            try:
                read = parse_sam_line(line)
                reference_end = calculate_reference_end(read)
                # print("Read name:", read.query_name)
                # print("Reference name:", read.reference_name)
                # print("Reference start position:", read.reference_start)
                # print("Reference end position:", reference_end)
                # print("CIGAR string:", read.cigarstring)
                # print("Sequence:", read.query_sequence)
                # print("Quality:", read.query_qualities)
                # print("Flags:", read.flag)
                # print("--------------")
                match_length=abs(int(reference_end)-int(read.reference_start))
                tmp_outs=[read.query_name,read.flag,read.reference_name,read.reference_start,read.mapping_quality,
                          read.cigarstring,'*',0,0,read.query_sequence,read.query_qualities,'*','*','*',reference_end,str(match_length)]
                tmp1.append(tmp_outs)
                # f1.write("\t".join([str(i) for i in tmp_outs])+"\n")
                # break
            except ValueError as e:
                print(f"Error parsing line: {e}", file=sys.stderr)
                continue

with open('test1.sam','w') as f1:
    for j in tmp1:
        f1.write('\t'.join([str(s) for s in j])+"\n")
