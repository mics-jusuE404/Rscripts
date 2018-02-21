## Because seqinr::write.fasta never gave me the result I wanted, borrowed this code from
## https://stackoverflow.com/questions/23374100/convert-table-into-fasta-in-r
## to write sequences and its name to a (multi)fasta file on disk:

Write.Fasta_AT <- function(Sequences, Names, FastaPath){
  if (length(Sequences) != length(Names)) stop="Sequences and Names of unequal length!"
  X <- data.frame(paste(">", Names, sep=""), Sequences)
  D <- do.call(rbind, lapply(seq(nrow(X)), function(i) t(X[i, ])))
  write.table(x = D, sep="\n", col.names = F, row.names = F, quote = F, file = FastaPath)
}

## Example:
Write.Fasta_AT(Sequences = my.seqs, Names = my.names, FastaPath = "/Path/To/file.fa")
