# Analysis decisions

## Deep intronic / splice variants (2026-10-05)
Precomputed genome-wide SpliceAI unavailable, and dbNSFP covers coding/splice-site SNVs only.
Three layers instead:
1. Wide candidate panel (~100 MVA / SAC / kinetochore / centromere genes): all rare intronic
   variants scored with SpliceAI and Pangolin; concordance between predictors reported.
2. "One hit" search: for any recessive gene with a single strong damaging variant and good
   phenotype fit, search introns, UTRs and SVs in that gene for a second allele.
3. Genome-wide sweep: rare (gnomAD AF < 0.1%) intronic variants in protein-coding genes,
   SpliceAI live on CPU (benchmark 1,000 variants first; GPU instance if > ~10 h),
   delta >= 0.2, intersected with Exomiser gene phenotype scores.
Optional: REMM in Exomiser for non-splice regulatory effects.
Limitation: splice effects are predictions; confirmation needs patient RNA or minigene assay.

## Submission construction (2026-10-06)
- Compound-het pairs built only from variants Exomiser scored under AR/XR modes.
- Hypervariable loci (HLA, KIR, MUC, olfactory receptors, immunoglobulin/TCR) excluded as
  mapping-artefact-prone in short-read WGS.
