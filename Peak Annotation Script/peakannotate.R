library(rtracklayer)
library(GenomicRanges)
library(dplyr)

# ----------------------------------------------------------------------------

# Import files
#annotation <- rtracklayer::import("v10DROPPED TriTrypDB-64_TbruceiLister427_2018.gff")
annotation <- readGFF("v10DROPPED TriTrypDB-64_TbruceiLister427_2018.gff")
tss <- readGFF("mapped_features_TSS.gff")

peaks <- rtracklayer::import("ChIP04_05_06_shared_peaks.bed")

#TSS <- annotation[annotation$type %in% c("dTSS", "sTSS")]
#TTS <- annotation[annotation$type %in% c("cTTS", "sTTS")]

# Convert to GRanges
annotation_gr <- as(annotation, "GRanges")
tss_gr <- as(tss,"GRanges")
peaks_gr <- as(peaks, "GRanges")
mcols(annotation_gr)$score <- NULL
annotation_cols <- names(mcols(annotation_gr))

# Ensure tss_gr has the same columns, adding missing ones as NA
tss_metadata <- mcols(tss_gr)

# Add NA-filled columns to tss_gr for those present in annotation_gr but not in tss_gr
for (col in setdiff(annotation_cols, names(tss_metadata))) {
  tss_metadata[[col]] <- NA
}

# Reorder columns to match annotation_gr
tss_metadata <- tss_metadata[, annotation_cols]

# Assign adjusted metadata back to tss_gr
mcols(tss_gr) <- tss_metadata

# Combine the GRanges objects
master_annotation_gr <- c(annotation_gr, tss_gr)

# Define the types you're interested in
higher_order_types <- c("protein_coding_gene", "pseudogene", "ncRNA_gene", "dTSS", "sTSS", "cTTS", "sTTS")

# Subset the annotation_gr object
high_order_annotation_gr <- master_annotation_gr[master_annotation_gr$type %in% higher_order_types]
print("prefilter:")
print(table(master_annotation_gr$type))
print("post filter:")
print(table(high_order_annotation_gr$type))



# Find all overlaps
#hits <- subsetByOverlaps(peaks_gr, annotation_gr)
hits <- findOverlaps(peaks_gr, high_order_annotation_gr)

print("Number of hits:")
print(length(unique(queryHits(hits))))

overlap_ranges <- pintersect(peaks_gr[queryHits(hits)], high_order_annotation_gr[subjectHits(hits)])
overlap_widths <- width(overlap_ranges)




# Create a data.frame to manipulate
overlap_df <- data.frame(
  peak_idx = queryHits(hits),
  feature_idx = subjectHits(hits),
  overlap_width = overlap_widths
)

best_hits_df <- overlap_df %>%
  group_by(peak_idx) %>%
  slice_max(order_by = overlap_width, n = 1, with_ties = FALSE) %>%
  ungroup()

#print(peaks_gr)
#print(best_hits_df)

# Subset peaks and annotations to unique best overlaps
unique_peaks <- peaks_gr[best_hits_df$peak_idx]
matched_annotations <- high_order_annotation_gr[best_hits_df$feature_idx]

# Attach annotation info to the peak metadata
mcols(unique_peaks) <- cbind(mcols(unique_peaks), mcols(matched_annotations))


# Export with 'type' and all other metadata included
# Rewrite the name column as feature type
# Set the name field in BED to something meaningful
mcols(unique_peaks)$name <- unique_peaks$type

# Identify the indices of peaks that *did* overlap
overlapped_idx <- unique(best_hits_df$peak_idx)

# Make a copy of the full peaks_gr
non_overlap_peaks <- peaks_gr

# Drop any rows whose index is in overlapped_idx
non_overlap_peaks <- non_overlap_peaks[-overlapped_idx]

# Ensure metadata columns exist and populate them
mcols(non_overlap_peaks)$type <- "intergenic"
mcols(non_overlap_peaks)$name <- "intergenic"

# Get the exact set of metadata columns from unique_peaks
common_cols <- names(mcols(unique_peaks))

# Add any missing columns to non_overlap_peaks, filling with NA
missing_cols <- setdiff(common_cols, names(mcols(non_overlap_peaks)))
for(col in missing_cols) {
  mcols(non_overlap_peaks)[[col]] <- NA
}

# Drop any extra columns (if there are any) and reorder to match unique_peaks
mcols(non_overlap_peaks) <- mcols(non_overlap_peaks)[, common_cols]

# Concatenate the GRanges
final_peaks <- c(unique_peaks, non_overlap_peaks)

# sort by genomic coordinate
final_peaks <- sort(final_peaks)
print(final_peaks)

# export your complete set
rtracklayer::export(final_peaks, "annotated_peaks.bed", format = "BED")