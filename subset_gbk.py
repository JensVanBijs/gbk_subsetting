#%% Imports

from Bio.SearchIO.BlastIO.blast_tab import BlastTabParser
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import re

SLICE_OFFSET = 10
BLASTDIR = "/home/jens/Desktop/6-fgtech-rasbora/utils/blast/output/hox"
SCAFFOLD = "./hoxda/hoxda-scaffold.fasta"

#%% Definitions

def blast_cluster_results(blast_dirpath:str, cluster:str) -> dict[str, tuple[list[int, int], int]]:
    """Process existing blast results into a machine-readable format.

    Args:
        blast_dirpath (str): Path to the directory containing all gene cluster BLAST results.
        cluster (str): Name of the directory containing the BLAST results for a single gene cluster.

    Returns:
        dict[str, tuple[list[int, int], int]]: Dictionary mapping gene name to gene location and strand.
    """
    collection_dict = dict()
    dirpath = os.path.join(blast_dirpath, cluster)
    genes = os.listdir(dirpath)
    ordered_genes = sorted(genes, key=lambda x: int(re.findall(r'\d+', x)[0]), reverse=False)

    for g in ordered_genes:
        blastpath = os.path.join(dirpath, g)
        parsed_blast = BlastTabParser(open(blastpath), comments=True)
        results = list(parsed_blast)[0]

        # HIT: get the collection of highest scoring hits from the blast results
        hit = list(results.__dict__["_items"].values())[0]

        # A High-scoring Segment Pair (HSP) is a local alignment with no gaps that achieves one of the highest alignment scores in a given search.
        hsp = hit[0]
        hit_range = hsp.hit_range
        strand = hsp.hit_strand

        collection_dict[g[:-4]] = (hit_range, strand)

    return collection_dict

def create_seq_obj(fasta_file:str, seq_id:str, name:str, desc:str) -> SeqIO.SeqRecord:
    """Create a documented Biopython sequence record object out of a fasta file.

    Args:
        fasta_file (str): Full filepath to a fasta file (containing a single sequence record).
        seq_id (str): Desired sequence ID.
        name (str): Desired sequence name.
        desc (str): Desired sequence description.

    Raises:
        AssertionError: FASTA file should contain exactly one sequence record

    Returns:
        SeqIO.SeqRecord: Biopython SeqRecord instance containing the sequence and all provided metadata.
    """
    parsed_seq = SeqIO.parse(fasta_file, "fasta")
    seqs = list(parsed_seq)

    if len(seqs) != 1:
        raise AssertionError("The FASTA file can only contain a single sequence: the contig/scaffold/chromosome containing the desired gene-cluster.")
        
    record = seqs[0]
    record.id = seq_id
    record.name = name
    record.desc = desc
    record.annotations = {"molecule_type": "DNA"}
    return record

def add_cds(seq_obj, name, start, end, strand):
    feature = SeqFeature(FeatureLocation(start, end), type="CDS", strand=strand, id=name)
    feature.qualifiers["gene"] = [name]
    # TODO: check if translation needs to be reversed if it is on the reverse strand
    feature.qualifiers["translation"] = [obj.seq[start:end].translate()]
    if strand == -1:
        feature.qualifiers["translation"] = feature.qualifiers["translation"][::-1]
    seq_obj.features.append(feature)
    return seq_obj

def subset_gbk(seq_obj, start, end):
    return seq_obj[start - SLICE_OFFSET : end + SLICE_OFFSET]

def save_gbk(seq_obj, filename):
    output_file = open(filename, "w")
    SeqIO.write(seq_obj, output_file, 'genbank')

#%%

blast_dirpath = os.path.abspath(BLASTDIR)
hoxda = blast_cluster_results(blast_dirpath, "hoxda")

obj = create_seq_obj(SCAFFOLD, "id123456", "hoxda_scaffold", "Scoffold containing the hoxda gene cluster from rasbora Lateristriata")
for key, value in hoxda.items():
    obj = add_cds(obj, key, value[0][0], value[0][-1], value[-1])

start = [x[0][0] for x in list(hoxda.values())]
end = [x[0][-1] for x in list(hoxda.values())]

save_gbk(subset_gbk(obj, min(start), max(end)), "./test.gbk")
