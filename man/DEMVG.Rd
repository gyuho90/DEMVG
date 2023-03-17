\name{DEMVG}
\alias{DEMVG}
\title{Calculates the DEMVG score}
\usage{
DEMVG(x,t_idx,c_idx)
}
\arguments{
\item{x}{DESeq normalized and vst transformed matrix}
\item{t_idx}{ The column index of treatment group}
\item{c_idx}{ The column index of control group }
}
\description{
DEMVG score measures differentially expressed genes (DEGs) in multivariate way.
}
\examples{

## do not RUN
# library(DESeq2)

# dds_bulkRNA <- DESeqDataSetFromMatrix( countData =  bulk_raw_count,
# colData = bulk_meta_data,  design =  ~ condition )

# vsd_bulkRNA <- DESeq2::varianceStabilizingTransformation(dds_bulkRNA)

# DEMVG_score<-DEMVG::DEMVG((assay(vsd_bulkRNA)),5:8,1:4)
}
