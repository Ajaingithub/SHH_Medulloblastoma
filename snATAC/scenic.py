## Running the GRN using scenic plus on the local machine
import os
import pandas as pd
os.chdir("/Users/abjain/Documents/UCSF/Projects/SHH/scATAc/GRN/")

out_dir = "outs"
os.makedirs(out_dir, exist_ok = True)

fragments_dict = {
    "10x_multiome_brain": "data/fragments.tsv.gz"
}

cell_data = pd.read_table("data/cell_data.tsv", index_col = 0)
cell_data.head()


chromsizes = pd.read_table(
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
chromsizes.head()

from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
os.makedirs(os.path.join(out_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)


bw_paths, bed_paths = export_pseudobulk(
    input_data = cell_data,
    variable = "VSN_cell_type",
    sample_id_col = "VSN_sample_id",
    chromsizes = chromsizes,
    bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    bigwig_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    path_to_fragments = fragments_dict,
    n_cpu = 10,
    normalize_bigwig = True,
    temp_dir = "/tmp",
    split_pattern = "-"
)

# We will need the paths to the bed files later on, so let's save them to disk.

with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv"), "wt") as f:
    for v in bw_paths:
        _ = f.write(f"{v}\t{bw_paths[v]}\n")

with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")

# Next we will use MACS to call peaks for each pseudobulk fragments.tsv.gz file.

bw_paths = {}
with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bw_paths.update({v: p})

from pycisTopic.pseudobulk_peak_calling import peak_calling
macs_path = "macs2"

os.makedirs(os.path.join(out_dir, "consensus_peak_calling/MACS"), exist_ok = True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.path.join(os.path.join(out_dir, "consensus_peak_calling/MACS")),
    genome_size = 'hs',
    n_cpu = 10,
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05,
    _temp_dir = '/tmp'
)



# Finally, it is time to derive the consensus peaks. To do so, we use the TGCA iterative peak filtering approach. First, each summit is extended a peak_half_width
# in each direction and then we iteratively filter out less significant peaks that overlap with a more significant one. During this procedure peaks will be merged
# and depending on the number of peaks included into them, different processes will happen:

#     1 peak: The original peak region will be kept
#     2 peaks: The original peak region with the highest score will be kept
#     3 or more peaks: The orignal peak region with the most significant score will be taken, and all the original peak regions in this merged peak region that
#  overlap with the significant peak region will be removed. The process is repeated with the next most significant peak (if it was not removed already) until 
#  all peaks are processed.

# This proccess will happen twice, first for each pseudobulk peaks; and after peak score normalization, to process all peaks together.


from pycisTopic.iterative_peak_calling import get_consensus_peaks
# Other param
peak_half_width=250
path_to_blacklist="/Users/abjain/Documents/UCSF/Resources/blacklist_region/hg38-blacklist.v2.bed"
# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = path_to_blacklist)

consensus_peaks.to_bed(
    path = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep =True,
    compression = 'infer',
    chain = False)

QC

# The next step is to perform QC for the scATAC-seq samples (in this case, only one run). There are several measurements and visualizations performed in this step:
    # Barcode rank plot
    # Duplication rate
    # Insertion size
    # TSS enrichment
    # Fraction of Reads In Peaks (FRIP)

# To calculate the TSS enrichment we need to provide TSS annotations. You can easily download them via the pycistopic tss get_tss command.
# In case you are unsure which column name is used by Ensembl to specify gene names in their databases, run the pycistopic tss gene_annotation_list and grep for your species.

# !mkdir -p outs/qc

# pycistopic tss get_tss \
#     --output outs/qc/tss.bed \
#     --name "hsapiens_gene_ensembl" \
#     --to-chrom-source ucsc \
#     --ucsc hg38

# Next, let’s calculate the QC metrics using the pycistopic qc command.
# pycistopic qc \
#     --fragments data/fragments.tsv.gz \
#     --regions outs/consensus_peak_calling/consensus_regions.bed \
#     --tss outs/qc/tss.bed \
#     --output outs/qc/10x_multiome_brain

# In case you have multiple samples, you can run the QC step in parallel as follows.

# regions_bed_filename = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
# tss_bed_filename = os.path.join(out_dir, "qc", "tss.bed")

# pycistopic_qc_commands_filename = "pycistopic_qc_commands.txt"

# # Create text file with all pycistopic qc command lines.
# with open(pycistopic_qc_commands_filename, "w") as fh:
#     for sample, fragment_filename in fragments_dict.items():
#         print(
#             "pycistopic qc",
#             f"--fragments {fragment_filename}",
#             f"--regions {regions_bed_filename}",
#             f"--tss {tss_bed_filename}",
#             f"--output {os.path.join(out_dir, "qc")}/{sample}",
#             sep=" ",
#             file=fh,
#         )

# Run the following command
# cat pycistopic_qc_commands.txt | parallel -j 4 {}

# Finally, we can visualize sample level statistics.
# These include:
# Barcode rank plot: The barcode rank plot shows the distribution of non-duplicate reads and which barcodes were inferred to be associated with cells. A steep drop-off 
# (‘knee’) is indicative of good separation between the cell-associated barcodes and the barcodes associated with empty partitions.
# Insertion size: ATAC-seq requires a proper pair of Tn5 transposase cutting events at the ends of DNA. In the nucleosome-free open chromatin regions, many molecules of Tn5
# can kick in and chop the DNA into small pieces; around nucleosome-occupied regions, and Tn5 can only access the linker regions. Therefore, in a good ATAC-seq library, you 
# should expect to see a sharp peak at the <100 bp region (open chromatin), and a peak at ~200bp region (mono-nucleosome), and other larger peaks (multi-nucleosomes). A 
# clear nucleosome pattern indicates a good quality of the experiment.

# Sample TSS enrichment: The TSS enrichment calculation is a signal to noise calculation. The reads around a reference set of TSSs are collected to form an aggregate 
# distribution of reads centered on the TSSs and extending to 1000 bp in either direction (for a total of 2000bp). This distribution is then normalized by taking the 
# average read depth in the 100 bps at each of the end flanks of the distribution (for a total of 200bp of averaged data) and calculating a fold change at each position 
# over that average read depth. This means that the flanks should start at 1, and if there is high read signal at transcription start sites (highly open regions of the genome)
# there should be an increase in signal up to a peak in the middle.

from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
import matplotlib.pyplot as plt

for sample_id in fragments_dict:
    fig = plot_sample_stats(
        sample_id = sample_id,
        pycistopic_qc_output_dir = "outs/qc"
    )

# The pycistopic qc command will determine automatic thresholds for the minimum number of unique number of fragments and the minumum TSS enrichment. In case you want to 
# change these thresholds or want to threhold based on FRIP, you can provide manually defined thresholds using the parameters: - unique_fragments_threshold - tss_enrichment_threshold
# - frip_threshold

# In this case we will use the automatically defined thresholds, please manualy inspect the quality metrics to make sure these thresholds are valid!    

# The barcode level statistics include:
# Total number of (unique) fragments
# TSS enrichment: The score at position in the TSS enrichmen score for for each barcode (at position 0, the TSS). Noisy cells will have a low TSS enrichment.
# FRIP: The fraction of reads in peaks for each barcode. Noisy cells have low FRIP values. However, this filter should be used with nuance, as it depends on the quality of the original peaks. For example, if there is a rare population in the sample, its specific peaks may be missed by peak calling algorithms, causing a decrease in their FRIP values.

from pycisTopic.qc import get_barcodes_passing_qc_for_sample
sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in fragments_dict:
    (
        sample_id_to_barcodes_passing_filters[sample_id],
        sample_id_to_thresholds[sample_id]
    ) = get_barcodes_passing_qc_for_sample(
            sample_id = sample_id,
            pycistopic_qc_output_dir = "outs/qc",
            unique_fragments_threshold = None, # use automatic thresholding
            tss_enrichment_threshold = None, # use automatic thresholding
            frip_threshold = 0,
            use_automatic_thresholds = True,
    )

for sample_id in fragments_dict:
    fig = plot_barcode_stats(
        sample_id = sample_id,
        pycistopic_qc_output_dir = "outs/qc",
        bc_passing_filters = sample_id_to_barcodes_passing_filters[sample_id],
        detailed_title = False,
        **sample_id_to_thresholds[sample_id]
    )

# Creating a cisTopic object
# In this step we will create a cisTopic object. This involves generating a count matrix containing fragment counts over consensus peaks (see above) for each cell barcode
# passing the QC metrices defined above. Blacklist regions will be removed from this count matrix.

path_to_regions = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
path_to_blacklist="/Users/abjain/Documents/UCSF/Resources/blacklist_region/hg38-blacklist.v2.bed"
pycistopic_qc_output_dir = "outs/qc"

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl

cistopic_obj_list = []
for sample_id in fragments_dict:
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[ sample_id_to_barcodes_passing_filters[sample_id] ]
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
        n_cpu = 1,
        project = sample_id,
        split_pattern = '-'
    )
    cistopic_obj_list.append(cistopic_obj)

# In this case we only have one sample, so only one cisTopic object has been generated. If you would have multiple samples, you would need to run the merge() function on your cisTopic object list.
cistopic_obj = cistopic_obj_list[0]
print(cistopic_obj)

import pickle
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)
