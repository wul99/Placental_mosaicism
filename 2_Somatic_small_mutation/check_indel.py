import pysam
import sys

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <bed> <bam1> <bam2> ... <output>")
        sys.exit(1)
    
    bed_file = sys.argv[1]
    bam_files = sys.argv[2:-1]  # bam file
    output = sys.argv[-1] # output
    
    # open file
    indel_data = []
    with open(bed_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                indel_data.append((chrom, pos, ref, alt))


    def process_bam(bam_file, chrom, start, ref, alt):
        has_mutation = "F"

        for read in bam_file.fetch(chrom, start-1, start + len(ref)):
            pos = read.pos
            for operation, length in read.cigartuples:
                if operation == 1:  # Insertion
                    if len(ref) < len(alt) and pos == start:
                        has_mutation = "T"
                        break
                    # if len(ref) < len(alt) and length == len(alt) - len(ref) and pos == start:
                    #     has_mutation = "T"
                    #     break
                    # if (len(ref) < len(alt) and read.query_sequence[start - read.pos - 1:start - read.pos + len(alt) - 1].upper() == alt and pos == start):
                    #     has_mutation = "T"
                    #     break
                elif operation == 2:  # Deletion
                    if len(ref) > len(alt) and pos == start:
                        has_mutation = "T"
                        break
                    pos += length  # Move reference position for deletions                
                    # if len(ref) > len(alt) and length == len(ref) - len(alt) and pos == start:
                    #     has_mutation = "T"
                    #     break
                    # pos += length  # Move reference position for deletions

                elif operation == 0 or operation == 7 or operation == 8:  # M, =, or X
                    pos += length  # Move reference position for matches/mismatches

                elif operation == 4:  # Soft clipping
                    pass

                elif operation == 5:  # Hard clipping
                    pass  # No need to update positions for hard clipping

                elif operation == 6:  # Padding
                    pass  # No need to update positions for padding

                elif operation == 3:  # Skipped region (N)
                    pos += length  # Move reference position for skipped region

        result=has_mutation
        return result

    results = []
    
    for chrom, start, ref, alt in indel_data:
        result_row = []

        for bam_path in bam_files:
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                result = process_bam(bam, chrom, start, ref, alt)
                result_row.append(str(result))
        
        result_line = [chrom, str(start), ref, alt] + result_row
        results.append("\t".join(result_line))

    # output
    with open(output, "w") as f:
        f.write("\n".join(results))


if __name__ == "__main__":
    main()
