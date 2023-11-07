import pandas as pd

##### JUST ADD TO pipeline2.py file in the place of bedtotrainset

def cistrometotrainset(filebed,filefasta,ENCODE_ID,TF):

    with open(filefasta) as f:  # fasta file to extract sequences
        fasta = f.readlines()
    seqs = [y.rstrip().upper().replace("N", "") for x, y in enumerate(fasta) if x % 2 != 0]
    peak_seqs = [x[(len(seqs[0]) // 2) - 30:(len(seqs[0]) // 2) + 30] for x in seqs]

    bed_data = pd.read_csv(filebed, sep="\t", header=None)  # bed narrowpeak file to extract scores
    bed_data.columns = ["chrom", "chromStart", "chromEnd", "peakindex", "score"]

    pbm_format = pd.DataFrame({0:bed_data["score"], 1:peak_seqs})
    pbm_format.sort_values(by=0,ascending=False,inplace=True)
    pbm_format.to_csv(f"outputs/ChIPseq_{ENCODE_ID}_{TF}.txt", header=None, index=False, sep="\t")
    return pbm_format
