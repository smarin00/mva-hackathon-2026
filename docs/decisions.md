# Analysis decisions

## Deep intronic / splice variants

**Plan (2026-10-05).** Precomputed genome-wide SpliceAI scores were unavailable
(the BaseSpace share is US-region only; the account is EU), and Exomiser's
dbNSFP scores cover coding and splice-site SNVs only. A three-layer contingency
was planned for the case that Exomiser found only one damaging allele in a
recessive gene: SpliceAI + Pangolin over rare intronic variants in the
candidate panel; a targeted search of introns, UTRs and SVs for a missing
second allele; and a genome-wide SpliceAI sweep of rare intronic variants.

**Outcome: not performed, because the contingency did not arise.** Exomiser
ranked two coding BUB1B alleles first as a compound-heterozygous pair
(p.Leu737\* and p.Asn1002Lys), confirmed by the Track 1 score. No second allele
had to be found.

**Targeted check instead (2026-10-19).** Because the Track 2 argument rests on
p.Asn1002Lys acting through protein stability, both alleles were scored with
SpliceAI 1.3.1 (GRCh38 annotation, +/-4,999 nt window) to exclude a splicing
mechanism (`make spliceai`). All delta scores <= 0.03 (threshold 0.2):

| Allele | AG | AL | DG | DL |
|---|---|---|---|---|
| p.Leu737\* (chr15:40209701 T>G) | 0.03 | 0.00 | 0.00 | 0.01 |
| p.Asn1002Lys (chr15:40220612 T>G), BUB1B | 0.02 | 0.00 | 0.00 | 0.00 |
| same position, PAK6 / read-through annotation | 0.02 | 0.00 | 0.01 | 0.00 |

No splice effect predicted for either allele: p.Leu737\* is null through NMD,
not exon skipping, and p.Asn1002Lys acts on the protein. Pangolin was not run.
Limitation: predictions only; RNA from the patient would confirm.

## Submission construction (2026-10-06)
- Compound-het pairs built only from variants Exomiser scored under AR/XR modes.
- Hypervariable loci (HLA, KIR, MUC, olfactory receptors, immunoglobulin/TCR) excluded as
  mapping-artefact-prone in short-read WGS.
