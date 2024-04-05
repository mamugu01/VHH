from Bio import SeqIO
from Bio import Align




def align_process(fastq_file, fasta_file):
    for query_id, query_seq in fasta_file.items():
        best_match = None
        best_score = 0

        # Iterate over each record in the FASTQ file
        for record in SeqIO.parse(fastq_file, "fastq"):
            # Perform pairwise alignment
            aligner = Align.PairwiseAligner()
            aligner.mode = "global"
            alignments = aligner.align(query_seq, record.seq)

            # Get the best alignment score
            top_alignment = alignments[0]
            alignment_score = top_alignment.score

            if alignment_score > best_score:
                best_score = alignment_score
                best_match = record

        # Output the best match
        if best_match is not None:
            print(f"Query ID: {query_id}")
            print(f"Best Match ID: {best_match.id}")
            print(f"Alignment Score: {best_score}")



def find_sequence(fastq_file, search_id):
    for record in SeqIO.parse(fastq_file, "fastq"):
        if record.id == search_id:
            print(f"Match found in {fastq_file}:")
            print(record.seq)
            return  

    # If no match is found
    print(f"No match found in {fastq_file}")

def main():
    fastq_file = "trim_1.fastqsanger"
    fasta_file = SeqIO.index("trim_DNA.fasta", "fasta")
    search_id = "M02058:630:000000000-LD8RJ:1:1111:10696:11653"
    align_process(fastq_file,fasta_file)
    #find_sequence(fastq_file, search_id)

main()
