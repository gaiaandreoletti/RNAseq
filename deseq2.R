#!/usr/bin/env Rscript

library(DESeq2)
library(tximport)
library(argparse)

parser = ArgumentParser()
parser$add_argument('--samplemap', help='Tab-separated sample map file', required=TRUE)
parser$add_argument('--countfiles_dir', help='directory containing countfiles', required=TRUE)
parser$add_argument('--tx2gene', help='transcript-to-gene file', required=FALSE)
parser$add_argument('--star_column', help='which column from the STAR ReadsPerGene.out.tab file to use', type='integer', default=4)
parser$add_argument('--out', help='output directory', default='.')
args = parser$parse_args()

sessionInfo()

load_data <- function(directory, tx2gene) {
    filenames = list.files(directory, recursive=TRUE)

    if ( all(grepl('quant.sf$', filenames)) ) {
        if (is.null(tx2gene)) {
            txOut = TRUE
        } else {
            txOut = FALSE
            sf = read.delim(paste0(directory, filenames[1]), header=TRUE, as.is=TRUE)
            print(paste("Length of quant.sf file:", nrow(sf)))
            print(paste("Length of tx2gene file:", nrow(tx2gene)))
            is_annotated_in_tx2gene = sf$Name %in% tx2gene[, 1]
            print(paste("Is annotated in tx2gene:"))
            print(table(is_annotated_in_tx2gene))
            if (!all(is_annotated_in_tx2gene)) {
                missing = sf$Name[!is_annotated_in_tx2gene]
                missing_df = data.frame(missing, missing)
                colnames(missing_df) = colnames(tx2gene)
                tx2gene = rbind(tx2gene, missing_df)
                print(paste("Added", length(missing), "transcripts to tx2gene."))
                print(missing)
            } else {
                print("All transcripts from quant.sf file are present in tx2gene file.")
            }
        }
        # do Salmon transcript import
        tx = tximport(paste0(directory, filenames), type='salmon', txOut=txOut, tx2gene=tx2gene, dropInfReps=TRUE)
        samplenames = sub('_quant.sf', '', filenames)
        samplenames = sub('[0-9]*/', '', samplenames)
        mat = as.matrix(tx$counts)
        colnames(mat) = samplenames
        # round to nearest integer
        mat = mat + 0.5
        mode(mat) = 'integer'
        return(mat)
    } else if ( all(grepl('quant.genes.sf$', filenames)) ) {
        # do Salmon transcript import
        tx = tximport(paste0(directory, filenames), type='salmon', txOut=TRUE, tx2gene=NULL, dropInfReps=TRUE)
        samplenames = sub('_quant.genes.sf', '', filenames)
        samplenames = sub('[0-9]*/', '', samplenames)
        mat = as.matrix(tx$counts)
        colnames(mat) = samplenames
        # round to nearest integer
        mat = mat + 0.5
        mode(mat) = 'integer'
        return(mat)
    } else if ( all(grepl('ReadsPerGene.out.tab.gz$', filenames))) {
        # do STAR gene import
        df.list = lapply(filenames, function (x) {
                temp = read.table(paste0(directory, '/', x), header = FALSE, row.names = 1, skip = 4)
                temp = temp[, args$star_column - 1, drop = FALSE]
                samplename = sub('.ReadsPerGene.out.tab.gz', '', x, fixed = TRUE)
                samplename = sub('[0-9]*/', '', samplename)
                colnames(temp) = samplename
                temp$gene = rownames(temp)
                return(temp)})
        df = Reduce(function(x, y) merge(x, y, by='gene', all.x=TRUE, all.y=TRUE), df.list)
        rownames(df) = df$gene
        df = df[, colnames(df) != 'gene']
        df[is.na(df)] = 0
        return(as.matrix(df))
    } else {
        system('dx-jobutil-report-error "countfiles not recognized as a supported type"')
    }
}

if ( ! is.null(args$tx2gene) ) {
    tx2gene <- read.delim(args$tx2gene, header=FALSE, as.is=TRUE)
    print("Loaded tx2gene.txt file")
    print(head(tx2gene))
} else {
    tx2gene <- NULL
}

countData <- load_data(directory=args$countfiles_dir, tx2gene)

sampleData <- read.table(args$samplemap, header = TRUE)

## Auto-magic factor ordering
## Any column that is a factor and has levels that have digits followed by two underscores
## (e.g. "01__cat") will have the levels automatically trimmed with the order of the levels preserved.
for (col in colnames(sampleData)) {
    if (is.factor(sampleData[, col]) & all(grepl('^[0-9]+__.+', levels(sampleData[, col])))) {
        print(paste('Factor level prefixes detected in column', col, '- stripping them now that order is set'))
        levels(sampleData[, col]) <- sub('^[0-9]+__(.+)', '\\1', levels(sampleData[, col]))
    }
}

sampleData <- sampleData[match(colnames(countData), sampleData$sample_id),]

dds = DESeqDataSetFromMatrix(countData, sampleData, design = ~ sample_id)
dds = estimateSizeFactors(dds)

save(dds, file=file.path(args$out, 'deseq2_object.RData'))

write.table(counts(dds, normalized=TRUE), file=file.path(args$out, 'deseq2_norm_count_matrix.txt'), sep = '\t', quote = FALSE)
