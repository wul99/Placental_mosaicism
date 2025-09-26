import argparse
import pysam
import os
import re  

parser = argparse.ArgumentParser(description="Somatic SNV phase.")
parser.add_argument("--sample_name", required=True)
parser.add_argument("--input_bam", "-b", help="bam file", required=True)
parser.add_argument("--somatic_vcf", "-s", help="unphased somatic vcf", required=True)
parser.add_argument("--germline_vcf", "-g", help="phased germline vcf (informative mutations)", required=True)
parser.add_argument("--output", "-o")

args = parser.parse_args()

sample = args.sample_name
input_bam = args.input_bam
somatic_vcf = args.somatic_vcf
germline_vcf = args.germline_vcf
output = args.output

def snv_indel_phase(sample_name,bamAlign_bam,somatic_vcf,germline_vcf,output):  
    somatic_vcf_read = pysam.VariantFile(somatic_vcf,"r")   
    for somatic_record in somatic_vcf_read:
        somatic_chrom = somatic_record.chrom
        somatic_pos = somatic_record.pos
        ref = somatic_record.ref
        alt = somatic_record.alts[0]

        query_names = []
        names = query_names.append

        bamAlign_f = pysam.AlignmentFile(bamAlign_bam,"rb")
        bamAlign_p = bamAlign_f.pileup(region=f"{somatic_chrom}:{somatic_pos}-{somatic_pos}",ignore_overlaps=False)
        for PileupColumn in bamAlign_p:
            if PileupColumn.reference_pos + 1 == int(somatic_pos):
                pileups = PileupColumn.pileups
                for PileupRead in pileups:
                    read = PileupRead.alignment
                    if not PileupRead.is_del and not PileupRead.is_refskip:
                        if len(alt) == len(ref) == 1:
                            if read.query_sequence[PileupRead.query_position].upper()==alt:
                                names(read.query_name)
                        else:
                            start = 0
                            for op, length in read.cigartuples:                                    
                                #insertion
                                if op == 1:
                                    insertion_base = read.query_sequence[start:start+length]
                                    if PileupRead.query_position + 1 == start and insertion_base.upper() == alt[1:] and len(alt)-len(ref) == length:
                                        names(read.query_name)
                                    start += length
                                #deletion
                                elif op == 2:
                                    if PileupRead.query_position + 1 == start and len(ref)-len(alt) == length:
                                        names(read.query_name)
                                elif op == 0 or op == 4 or op == 7 or op == 8:
                                    start += length
                                else:
                                    pass
        query_name = set(query_names)

        bamAlign_f2 = bamAlign_f.fetch(region=''.join([somatic_chrom,':',str(int(somatic_pos)-1000 if int(somatic_pos)-1000>0 else 1),'-',str(int(somatic_pos)+1000)]))
        tem_bam = "tem.bam"
        tem_bam_write = pysam.AlignmentFile(tem_bam, 'wb', template=bamAlign_f)
        for read in bamAlign_f2:
            if read.query_name in query_name:
                tem_bam_write.write(read)
        tem_bam_write.close()
        bamAlign_f.close()
        pysam.index(tem_bam)

        tem_bam_read = pysam.AlignmentFile(tem_bam, 'rb')
        F = []
        M = []         
        germline_vcf_read = pysam.VariantFile(germline_vcf,"r")
        for germline_record in germline_vcf_read:
            if len(germline_record.alts) == 1:
                germline_GT = germline_record.samples[sample_name].values()
                germline_chrom = germline_record.chrom
                germline_pos = int(germline_record.pos)
                if germline_GT[0][0] != germline_GT[0][1]:
                    if germline_GT[0][0] == 0:
                        Base_F = germline_record.ref
                        Base_M = germline_record.alts[0]
                    elif germline_GT[0][0] == 1:
                        Base_F = germline_record.alts[0]
                        Base_M = germline_record.ref
                    bamAlign_p2 = tem_bam_read.pileup(region=f"{germline_chrom}:{germline_pos}-{germline_pos}",ignore_overlaps=False)
                    for germline_PileupColumn in bamAlign_p2:
                        if germline_PileupColumn.reference_pos + 1 == germline_pos:
                            germline_pileups = germline_PileupColumn.pileups
                            for germline_PileupRead in germline_pileups:
                                germline_read = germline_PileupRead.alignment
                                if not germline_PileupRead.is_del and not germline_PileupRead.is_refskip:
                                    # snv                 
                                    if len(Base_F) == len(Base_M) == 1:
                                        if germline_read.query_sequence[germline_PileupRead.query_position].upper()==Base_F:
                                            F.append(germline_pos)
                                        elif germline_read.query_sequence[germline_PileupRead.query_position].upper()==Base_M:
                                            M.append(germline_pos)
                                    # indel
                                    else:
                                        start = 0
                                        for op, length in germline_read.cigartuples:
                                            # insertion
                                            if op == 1:
                                                insertion_base = germline_read.query_sequence[start:start+length]
                                                if PileupRead.query_position + 1 == start and len(germline_record.alts[0]) - len(germline_record.ref) == length: 
                                                    if insertion_base.upper() == Base_F[1:]:
                                                        F.append(germline_pos)
                                                    elif insertion_base.upper() == Base_M[1:]:
                                                        M.append(germline_pos)
                                                start += length
                                            # deletion
                                            elif op == 2:
                                                if PileupRead.query_position + 1 == start and len(germline_record.ref) - len(germline_record.alts[0]) == length:
                                                    if len(Base_F) == 1:
                                                        F.append(germline_pos)
                                                    elif len(Base_M) == 1:
                                                        M.append(germline_pos)
                                            elif op == 0 or op == 4 or op == 7 or op == 8:
                                                start += length
                                            else:
                                                pass
        with open(output, 'a') as phase_txt_write:
            if len(F) != len(M):
                if len(F) == 0:
                    phase_txt_write.write(line.rstrip('\n')+"\tM\n")
                elif len(M) == 0:
                    phase_txt_write.write(line.rstrip('\n')+"\tF\n")
                else:
                    phase_txt_write.write(line.rstrip('\n')+"\tFM\n")
                    print(sample_name+"\t"+line.rstrip('\n')+"\tFM\t"+str(F)+"|"+str(M)+"\n")
            elif len(F) == len(M):
                if len(F) > 0:
                    phase_txt_write.write(line.rstrip('\n')+"\tFM\n")
                    print(sample_name+"\t"+line.rstrip('\n')+"\tFM\t"+str(F)+"|"+str(M)+"\n")
                elif len(F) == 0:
                    phase_txt_write.write(line.rstrip('\n')+"\tunphase\n")
        germline_vcf_read.close()
        tem_bam_read.close()
        try:  
            if os.path.exists(tem_bam):  
                os.remove(f"{tem_bam}")
                os.remove(f"{tem_bam}.bai")
        except OSError as e:  
            print(f"Error occurred while deleting temporary file {tem_bam}: {e.strerror}")

snv_indel_phase(sample_name,input_bam,somatic_vcf,germline_vcf,output)


