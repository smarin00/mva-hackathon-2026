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

## Track 1 rules generalised (October 2026)
- Lessons from the gene-level review became rules applied to every genome
  (`bin/make_submission.py`): population frequency must fit the mode of
  inheritance (at most 0.1% for AD/XD, 2% for AR/XR, cf. ACMG BS1), and genes
  with no OMIM or Orphanet disease in the HPO annotations (v2026-09-01) are
  ranked last rather than removed. On this genome the rules alone produce the
  same four submitted variants; the review only moves LZTR1 to a secondary
  finding (`config/cases/mva1_proband.curation.tsv`).

## Final Track 1 submission (2026-10-24)
- The read-level submission (`findings/track1/submission_full.csv`) replaces the
  first one (provided VCF, `findings/track1/submission_vcf.csv`). Both contain the same
  four variants; the read-level route adds independent calling, phasing and the
  mosaic, structural, repeat and aneuploidy checks.

## Reproducibility check (2026-10-22)
- A clean clone was run end to end with `make all` (about 15 h). The submission
  and every Track 2 table were identical to the committed ones. Aneuploidy
  estimates differed by at most 0.03 percentage points: the mosdepth coverage
  differed, because bwa-mem without a fixed batch size (`-K`) batches reads by
  thread timing. Adding `-K` would make alignment deterministic, but was not
  done before the deadline because it could not be re-validated on the full
  run; it is listed as a limitation.

## Track 2 statistics and counts (2026-10-22)
- ClinVar check of the atlas: the expected share of each bin is taken from
  missense changes reachable by a single nucleotide substitution (the 6,173
  changes AlphaMissense scores), not from all 19 substitutions per residue,
  because ClinVar can only hold the former. Benign variants in the damaging bin:
  1 of 33 against 34% expected (one-sided binomial p = 2e-5; tolerated bin 19 of
  33 against 26%, p = 1e-4). The earlier all-substitution baseline (43%,
  p = 2e-7) overstated the result.
- No AUROC is reported for the ClinVar check: ClinVar holds a single pathogenic
  BUB1B missense outside the calibration set.
- Drug counts in the screen are parent compounds (salt forms merged): 15 in the
  desired direction (12 HDAC inhibitors, 2 BCL-xL inhibitors, PRI-724; all
  curated out) and 31 do-no-harm (Aurora B 18, HSP90 9 including CD24Fc,
  NAMPT 2, MPS1 2).
- The Aurora B entry of the do-no-harm list is marked contested: inhibiting
  Aurora B restored chromosome congression in cells whose BubR1 cannot recruit
  PP2A-B56 (Xu et al., Biol Open 2013). It stays on the list at oncology doses;
  low-dose inhibition is an exploratory step of the validation plan.
- The child presented with rhabdomyosarcoma (HP:0002859). The cancer-treatment
  check is therefore framed around the child's own treatment, and two case
  reports of poor chemotherapy tolerance in MVA (Ochiai 2019; Nishitani-Isa
  2019) are cited.
