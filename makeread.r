library(polyester)
library(Biostrings)

# read the quant.sf
makeReads <- function(quant_file, fastapath, output_dir, repeat_times=1){

    quants <- read.table(quant_file, header = T, stringsAsFactors = F)

    fasta <- readDNAStringSet(fastapath)

    if (!all(quants$Transcript == names(fasta))) {
        stop("quant assignment Not match with fasta file ")
    }
    NumRead <- quants$NumRead

    readscount <- matrix(rep(NumRead, repeat_times)
                         ,nrow = length(NumRead)
                         ,ncol = repeat_times)
    
    simulate_experiment_countmat(fasta = fastapath, 
                                 readmat = readscount,
                                 outdir = output_dir)
    
}


makeReads(
    quant_file = "mixed.txt",
          fastapath = "mixed.fa",
    output_dir = "estimatedReadRaw",
    repeat_times = 5
)

sink()
q(save = "no")

