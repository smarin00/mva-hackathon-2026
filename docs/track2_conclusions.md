# Track 2: Conclusions so far

Working record of the mechanistic findings that Track 2 builds on. Each point
is tagged **[data]** (from this child's sequencing), **[literature]**
(published experimental evidence), or **[hypothesis]** (our inference, to be
tested or stated as such).

Status: 19 Oct 2026. Pending: ESM-2 calibrated score for N1002K.

---

## 1. Genotype

| Allele | Position (GRCh38) | cDNA | Protein | Consequence |
|---|---|---|---|---|
| 1 | chr15:40209701 T>G | c.2210T>G | p.Leu737\* | stop-gain |
| 2 | chr15:40220612 T>G | c.3006T>G | p.Asn1002Lys | missense, pseudokinase domain |

Transcript: ENST00000287598 (MANE Select), 23 exons, 1049 aa.

- **[data]** Both variants were called from reads by DeepVariant and in the
  vendor VCF; Exomiser ranks BUB1B first (gene score 0.587 / 0.554).
- **[data]** Phase is **unresolved**. The variants are 10.9 kb apart, beyond
  the reach of short-read pairs, and WhatsHap found no shared phase block.
  Being in trans is required for the diagnosis: in cis, the other copy would
  be wild type and the child a carrier.
- **[literature]** The evidence for trans is phenotypic: MVA with
  constitutional aneuploidy arises only from biallelic BUB1B loss (Hanks
  2004). Parental testing or long reads would settle it; neither is
  available here.

## 2. L737\* is a null allele

- **[data]** NMD geometry (`bin/nmd_check.py`): the PTC lies in exon 17 of
  23, **747 nt upstream of the last exon–exon junction**, about 15 times the
  50-nt threshold. Verdict: NMD-sensitive.
- **[literature]** Truncating MVA alleles are not detectable as protein in
  patient cells, and their transcripts are reduced or absent (Suijkerbuijk
  2010).
- Note: a clean truncation at 731 fully rescues checkpoint function
  (Suijkerbuijk 2012). The 737\* product *would* be functional if made. NMD
  is what makes this allele null, not the loss of protein sequence.

**Conclusion: all of this child's BubR1 comes from the N1002K allele.**

## 3. The disease mechanism is protein dosage

- **[literature]** No MVA patient carries two truncating alleles. Biallelic
  null is likely embryonic lethal, so MVA1 requires residual BubR1 (Hanks
  2004).
- **[literature]** The BubR1 kinase domain is a **pseudokinase**:
  catalytically inactive and dispensable for the spindle checkpoint
  (Suijkerbuijk 2012).
- **[literature]** MVA missense alleles in that domain (I909T, L1012P,
  R727C) are **stability mutants**: about 2-fold faster turnover, dependent
  on HSP90, cleared by the proteasome, 2–6-fold lower abundance in patient
  cells. Restoring their level restores function (Suijkerbuijk 2010).
- **[literature]** The dose response is steep. About 6% of normal BubR1
  causes missegregation in most cells; about 13% has virtually no effect
  (Suijkerbuijk 2010). Mice with ~10% BubR1 survive with aneuploidy and
  growth defects.

## 4. N1002K: a candidate stability mutant

- **[data]** N1002K lies in the pseudokinase domain, **10 residues from
  L1012P**, an established HSP90-dependent stability mutant.
- **[hypothesis]** N1002K destabilises BubR1 by the same mechanism and
  leaves the child near the ~6–13% threshold.
- **Test (running):** ESM-2 masked-marginal score, calibrated against six
  MVA alleles with known experimental stability classes
  (`bin/esm_variant_score.py`). If the model separates the classes, N1002K's
  placement is evidence; if not, it is supporting context only.

## 5. Therapeutic implications

- **[hypothesis]** This is a proteostasis problem with a narrow window:
  raising BubR1 roughly 2-fold, from near 6% toward 13%, may be enough. The
  aim is to stabilise the N1002K protein, not to restore enzyme activity
  (there is none).
- **[hypothesis] Do-no-harm:** **HSP90 inhibitors** would be expected to
  collapse the remaining BubR1. Mutant BubR1 is an HSP90 client, and HSP90
  inhibition depleted it in patient-derived cells (Suijkerbuijk 2010). Flag
  for any future oncology treatment in a cancer-predisposed child.
- **[literature]** Proteasome inhibition rescued mutant BubR1 levels in
  vitro. That is a proof of mechanism, not a treatment suggestion.
- Not relevant: read-through or NMD-modulating strategies for L737\*. They
  are theoretically interesting because the truncated product would be
  functional, but they are a far weaker lever than stabilising the allele
  that already makes protein.

## 6. Aneuploidy from bulk sequencing depth: unconfirmed

- **[data]** GC-corrected depth shows chr16, 17 and 19 elevated by 1.3–1.9%.
  The shift persists under three GC-correction models, but chromosome GC
  content predicts the depth ratio under all three (Spearman ρ 0.48–0.69).
  BAF shows no corresponding signal.
- **[data]** These are reported as **flagged, unconfirmed**, consistent with,
  but not proven to be, residual GC/replication-timing bias. Practical
  detection limit: above ~3% mosaicism.
- **[literature]** This is expected. MVA aneuploidy is random across cells
  and averages out in bulk DNA; the cytogenetic diagnosis rests on the
  clinical karyotype.

| GC model | ρ | chr16 | chr17 | chr19 |
|---|---|---|---|---|
| local | +0.65 | +1.30% | +1.51% | +1.92% |
| 2D, leave-one-out | +0.48 | +1.25% | +1.42% | +1.91% |
| poly, leave-one-out | +0.69 | +1.46% | +1.69% | +2.28% |

## 7. Open items

- ESM-2 calibrated score for N1002K
- Structural check of N1002K in the AlphaFold model (buried vs. surface,
  local contacts)
- Candidate generation: proteostasis modulators, drug–target evidence
  (Open Targets), do-no-harm list beyond HSP90
- Pitch video (3 min)

## References

- Hanks S. et al. Constitutional aneuploidy and cancer predisposition caused
  by biallelic mutations in BUB1B. *Nat Genet* 36, 1159–1161 (2004).
- Suijkerbuijk S.J.E. et al. Molecular causes for BUBR1 dysfunction in the
  human cancer predisposition syndrome mosaic variegated aneuploidy.
  *Cancer Res* 70, 4891–4900 (2010).
- Suijkerbuijk S.J.E. et al. The vertebrate mitotic checkpoint protein BUBR1
  is an unusual pseudokinase. *Dev Cell* 22, 1321–1329 (2012).
- BubR1 TPR domain supports the mitotic checkpoint by promoting MCC formation
  and MCC-APC/C interaction. *Commun Biol* (2026),
  doi:10.1038/s42003-026-10384-9.
