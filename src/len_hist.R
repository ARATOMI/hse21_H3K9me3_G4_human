
source('lib.R')

###

#NAME <- 'G4_Li_KPDS'
#NAME <- 'H3K9me3_A549.ENCFF164FDB.hg19'
#NAME <- 'H3K9me3_A549.ENCFF164FDB.hg38'
#NAME <- 'H3K9me3_A549.ENCFF494QKI.hg19'
#NAME <- 'H3K9me3_A549.ENCFF494QKI.hg38'

###

bed_df <- read.delim(paste0(DATA_DIR, NAME, '.bed'), as.is = TRUE, header = FALSE)
#colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
colnames(bed_df) <- c('chrom', 'start', 'end')
bed_df$len <- bed_df$end - bed_df$start

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('len_hist.', NAME, '.pdf'), path = OUT_DIR)

