import pysam
import sys
from collections import Counter

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <input> <bam1> <bam2> ... <output>")
        sys.exit(1)
    
    input = sys.argv[1] #input
    bam_files = sys.argv[2:-1]  #bam file
    output = sys.argv[-1] #output
    
    sv_data = []

    # open file
    with open(input, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) == 6:
                chrom = parts[0]
                start1 = int(parts[1])
                end1 = int(parts[2])
                start2 = int(parts[4])
                end2 = int(parts[5])
                sv_data.append((chrom, start1, end1, start2, end2))


    def find_breakpoint(bam_file, chrom, bp1, bp2):
        pos_list = []
        if bp1 <= bp2:
            p1=bp1
            p2=bp2
        else:
            p1=bp2
            p2=bp1
        for read in bam_file.fetch(chrom, p1-10, p2+10):
            pos = read.pos
            for operation, length in read.cigartuples:
                if operation == 1:  # Insertion
                    pass
                elif operation == 2:  # Deletion
                    pos += length  # Move reference position for deletions

                elif operation == 0 or operation == 7 or operation == 8:  # M, =, or X
                    pos += length  # Move reference position for matches/mismatches
                    
                elif operation == 4:  # Soft clipping
                    pos_list.append(pos)

                elif operation == 5:  # Hard clipping
                    pass  # No need to update positions for hard clipping

                elif operation == 6:  # Padding
                    pass  # No need to update positions for padding

                elif operation == 3:  # Skipped region (N)
                    pos += length  # Move reference position for skipped region

        pos_counter = Counter(pos_list)

        # find the position with the highest frequency
        if pos_counter: 
            most_common_pos = pos_counter.most_common(1)[0][0]
            if pos_counter.most_common(1)[0][1] >= 2 :
                return most_common_pos


    results = []


    for chrom, start1, end1, start2, end2 in sv_data:
        all_sv = []
        for bam_file_path in bam_files:
            bam_file = pysam.AlignmentFile(bam_file_path, "rb")
            left_bp = find_breakpoint(bam_file, chrom, start1, start2)
            right_bp = find_breakpoint(bam_file, chrom, end1, end2)
            if left_bp is not None and right_bp is not None:
                sv=chrom+":"+str(left_bp)+"-"+str(right_bp)
            else:
                sv="null"
            all_sv.append(str(sv))
            bam_file.close()

        pos_counter = Counter(pos for pos in all_sv if pos != "null")
        if pos_counter and (len(pos_counter)==1 or (len(pos_counter)>1 and pos_counter.most_common(1)[0][1] > 1)):
            result = pos_counter.most_common(1)[0][0]
        else:
            result = chrom+":"+str(start1)+"-"+str(end1)+"\t"+chrom+":"+str(start2)+"-"+str(end2)

        results.append((result))

    # output
    with open(output, "w") as f:
        for item in results:
            f.write(item + '\n')

if __name__ == "__main__":
    main()
