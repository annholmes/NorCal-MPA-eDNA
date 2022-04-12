#install qiime2
wget
conda install wget
wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-osx-conda.yml
conda env create -n qiime2-2022.2 --file qiime2-2022.2-py38-osx-conda.yml
rm qiime2-2022.2-py38-osx-conda.yml
conda activate qiime2-2022.2
qiime --help #test installation

#create MultiplexedPairedEndBarcodeInSequence qza file using muxed-pe-barcode-in-seq folder
#mux-pe-barcode-in-seq folder contains forward and reverse sequence fastq.gz files 
mkdir muxed-pe-barcode-in-seq
qiime tools import \
  --type MultiplexedPairedEndBarcodeInSequence \
  --input-path muxed-pe-barcode-in-seq \
  --output-path multiplexed-seqs.qza
  
#demultiplex with cutadapt
# --i-seqs The paired-end sequences to be demultiplexed
# --m-forward-barcodes-column sample metadata column listing the per-sample barcodes for the forward reads
# --m-reverse-barcodes-column sample metadata column listing the per-sample barcodes for the reverse reads   
# --o-per-sample-sequences demux-full.qza The resulting demultiplexed sequences
# --o-untrimmed-sequences unmatched.qza The sequences that were unmatched to barcodes

qiime cutadapt demux-paired \
  --i-seqs multiplexed-seqs.qza \
  --m-forward-barcodes-file sample-metadata.tsv \
  --m-forward-barcodes-column forward-barcodes \
  --m-reverse-barcodes-file sample-metadata.tsv \
  --m-reverse-barcodes-column reverse-barcodes \
  --o-per-sample-sequences demux-full.qza \
  --o-untrimmed-sequences unmatched.qza

# Summarize/visualize demultiplexed output 
# File can be viewed on https://view.qiime2.org/
qiime demux summarize \
  --i-data demux-full.qza \
  --o-visualization demux-full.qza.qzv
  
# export demultiplexed files as fastq.gz files 
mkdir extracted-demux-full
qiime tools extract \
  --input-path demux-full.qza \
  --output-path extracted-demux-full

