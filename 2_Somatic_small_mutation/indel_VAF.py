import pysam
import sys

def main():
    if len(sys.argv) < 4:
        print("Usage: python script.py <bed> <bam1> <bam2> ... <output_vaf> <output_depth>")
        sys.exit(1)
    
    bed_file = sys.argv[1]
    bam_files = sys.argv[2:-2]  # bam file
    output_vaf = sys.argv[-2]   # output_vaf
    output_depth = sys.argv[-1] # output_depth
    
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
        count = 0
        depth = 0
        
        for read in bam_file.fetch(chrom, start-1, start):
            depth += 1
            pos = read.pos
            for operation, length in read.cigartuples:
                if operation == 1:  # Insertion
                    if len(ref) < len(alt) and length == len(alt) - len(ref) and pos == start:
                        count += 1
                    elif (len(ref) < len(alt) and read.query_sequence[start - read.pos - 1:start - read.pos + len(alt) - 1].upper() == alt and pos == start):
                        count += 1
                
                elif operation == 2:  # Deletion
                    if len(ref) > len(alt) and length == len(ref) - len(alt) and pos == start:
                        count += 1
                    pos += length
                
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
        
        vaf = 0.0 if depth == 0 else count / depth
        return vaf, depth


    vaf_results = []
    depth_results = []
    
    for chrom, start, ref, alt in indel_data:
        vaf_row = []
        depth_row = []
        key = f"{chrom}_{start}_{ref}_{alt}"
        
        for bam_path in bam_files:
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                vaf, depth = process_bam(bam, chrom, start, ref, alt)
                vaf_row.append(str(vaf))
                depth_row.append(str(depth))
        
        vaf_results.append("\t".join([key] + vaf_row))
        depth_results.append("\t".join([key] + depth_row))

    # output
    with open(output_vaf, "w") as f:
        f.write("\n".join(vaf_results))
    
    with open(output_depth, "w") as f:
        f.write("\n".join(depth_results))

if __name__ == "__main__":
    main()
