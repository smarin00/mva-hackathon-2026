# Track 2: Conclusions so far

Working record of the mechanistic findings that Track 2 builds on. Each point
is tagged **[data]** (from this child's sequencing), **[literature]**
(published experimental evidence), or **[hypothesis]** (our inference, to be
tested or stated as such).

Status: 24 Oct 2026. Final; the submitted report is
`docs/smarin00_track2_report.pdf`, candidates in `track2_candidates.md`.

The child presented with rhabdomyosarcoma (HP:0002859), the event that triggered
the genetic investigation.

---

## 1. Genotype

| Allele | Position (GRCh38) | cDNA | Protein | Consequence |
|---|---|---|---|---|
| 1 | chr15:40209701 T>G | c.2210T>G | p.Leu737\* | stop-gain |
| 2 | chr15:40220612 T>G | c.3006T>G | p.Asn1002Lys | missense, pseudokinase domain |

Transcript: ENST00000287598 (MANE Select), 23 exons, 1050 aa.

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

- **[data]** SpliceAI predicts no splice effect for either allele (all delta
  scores <= 0.03; `make spliceai`, `docs/decisions.md`). L737\* is null through
  NMD, not exon skipping; N1002K acts on the protein, not on splicing.

**Conclusion: all of this child's BubR1 comes from the N1002K allele.**

## 3. The disease mechanism is protein dosage

- **[literature]** No MVA patient carries two truncating alleles. Biallelic
  null is likely embryonic lethal, so MVA1 requires residual BubR1 (Hanks
  2004).
- **[literature]** The BubR1 kinase domain is a **pseudokinase**:
  catalytically inactive and dispensable for the spindle checkpoint
  (Suijkerbuijk 2012).
- **[literature]** MVA missense alleles in or packed against that domain
  (I909T, L1012P, R727C; R727 lies just before the curated boundary at 766 but
  contacts Trp803 and Phe912 in the AlphaFold model) are **stability mutants**: about 2-fold faster turnover, dependent
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
- **[data]** ESM-2 (650M) masked-marginal score, calibrated against six MVA
  alleles with experimentally measured stability classes
  (`bin/esm_variant_score.py`):

  | Variant | Experimental class | ESM LLR | Percentile in domain |
  |---|---|---|---|
  | L1012P | destabilising | −6.94 | 12.8 |
  | **N1002K** | **patient** | **−5.81** | **18.2** |
  | I909T | destabilising | −5.71 | 18.8 |
  | L844F | destabilising | −3.33 | 34.8 |
  | R727C | destabilising | −2.80 | 40.4 |
  | Q921H | stable | −0.24 | 78.7 |
  | R550Q | stable | +0.37 | 87.3 |

  The model separates the classes completely, and N1002K falls between
  L1012P and I909T, the two alleles for which HSP90 dependence and rescue by
  overexpression were shown directly.
- Caveats: n = 6 (a perfect 4-vs-2 split occurs by chance with p ≈ 0.07);
  ESM scores evolutionary constraint, not folding, and is a stability proxy
  here only because the pseudokinase domain's known role is structural;
  R550Q lies outside the domain, so Q921H is the only in-domain stable
  control.
- **[data]** Structural context in the AlphaFold model (AF-O60566-F1;
  `bin/structure_context.py`). N1002 is confidently modelled (pLDDT 91) and
  partially buried (relative solvent accessibility 0.16) in an aromatic pocket
  (Trp973, Phe977, Trp978, Val998). Its side chain makes three hydrogen bonds
  to backbone atoms of Trp978 (N, O) and Val998 (O), bridging two segments of
  the fold. Lysine is longer and charged: it cannot reproduce this
  donor-acceptor bridge and places a positive charge in an aromatic pocket.
- **[data] Calibration of the structure:** burial does **not** discriminate the
  experimental classes. The verified destabilisers are buried (RSA 0.00–0.06),
  but so is the stable control Q921H (0.01), a conservative polar swap; R550Q
  is in a disordered region (pLDDT 35) and cannot be assessed. The structure
  therefore supplies a physical mechanism, not independent discrimination.
  N1002 and L1012 do not share contacts, so we claim the same domain and
  sequence region, not the same pocket.
- Monomer model: interfaces with Bub3, Mad2 and Cdc20 are not represented.

| Line of evidence | Result | Weight |
|---|---|---|
| ESM-2, calibrated | Classes separate; N1002K between L1012P and I909T | Moderate (n = 6, p ≈ 0.07) |
| Structure | Confident, partially buried; breaks a side-chain H-bond bridge | Physical mechanism; burial alone does not discriminate |
| AlphaMissense | N1002K 0.92, likely pathogenic; does not separate the experimental stability classes (R727C, L844F called likely benign) | Independent convergence on the patient allele; weaker than calibrated ESM-2 for mild stability alleles |
| Literature | Same domain as four verified stability mutants; domain role is structural | Strong prior |

- **Conclusion [hypothesis, supported by four lines of evidence]:** N1002K
  most likely destabilises BubR1. No single line is proof; experimental
  confirmation is the cycloheximide-chase step of the validation plan in
  `track2_candidates.md`.

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
- **[literature]** BubR1 acetylation at K250 is mitosis-only and protects it
  from APC/C degradation until checkpoint silencing, when HDAC2/3 remove it
  (Choi 2009; Park 2017). It sets mitotic timing, not steady-state level
  (Park 2013). HDAC inhibition (valproate) is therefore an exploratory
  hypothesis, not a dose-raising candidate; see `track2_candidates.md`.
- Not relevant: read-through or NMD-modulating strategies for L737\*. They
  are theoretically interesting because the truncated product would be
  functional, but they are a far weaker lever than stabilising the allele
  that already makes protein.

## 6. Aneuploidy from bulk sequencing depth: unconfirmed

- **[data]** GC-corrected depth shows chr16, 17 and 19 elevated by 1.3–2.3%,
  depending on the GC model (table below).
  The shift persists under three GC-correction models, but chromosome GC
  content predicts the depth ratio under all three (Spearman ρ 0.48–0.69).
  BAF shows no corresponding signal.
- **[data]** These are reported as **flagged, unconfirmed**, consistent with,
  but not proven to be, residual GC/replication-timing bias. Noise-based
  detection limit: about 2% of cells (`findings/aneuploidy/*.limit`).
- **[literature]** This is expected. MVA aneuploidy is random across cells
  and averages out in bulk DNA; the cytogenetic diagnosis rests on the
  clinical karyotype.

| GC model | ρ | chr16 | chr17 | chr19 |
|---|---|---|---|---|
| local | +0.65 | +1.30% | +1.51% | +1.92% |
| 2D, leave-one-out | +0.48 | +1.25% | +1.42% | +1.91% |
| poly, leave-one-out | +0.69 | +1.46% | +1.69% | +2.28% |

## 7. BUB1B variant atlas and ClinVar check

- **[data]** The atlas (`findings/track2/atlas/`) recovers the functional map
  from sequence alone: KEN1 (mean ESM LLR −11.2, more constrained than 98.6% of
  positions), KEN2 (−9.0), GLEBS (−7.8) and TPR (−6.9) are the most constrained
  features; the disordered linker 450–650 (median pLDDT 34) is nearly
  unconstrained (−1.0).
- **[data]** Premature stops at codons 1–969 are NMD-sensitive; codons
  970–1050 escape NMD and keep all checkpoint-essential features.
- **[data]** ClinVar (`findings/track2/clinvar_validation.tsv`): 1,494 BUB1B
  missense records, 97% VUS/conflicting. Benign variants (33 unique) fall in
  the damaging bin at 3% against 34% of single-nucleotide-reachable changes
  (one-sided binomial p = 2e-5), and in the tolerated bin at 58% against 26%
  (p = 1e-4). (An earlier baseline over all 19 substitutions, 43% and 20%,
  overstated this; see `docs/decisions.md`.) One non-calibration pathogenic
  variant (Q467H) only, so the pathogenic side cannot be tested.
- **[data]** Of 1,427 unique VUS, the atlas places 107 in the stability-defect
  class and 419 as likely tolerated; AlphaMissense and ESM-2 agree on 12%
  (damaging) and 28% (tolerated) of VUS records.
- **[data]** ClinVar lists the verified destabilisers L1012P, L844F and R727C,
  and the patient's N1002K, as VUS.

## 8. Related genes

- **[data]** From the HPO annotations the pipeline finds six other genes linked
  to MVA: CEP57 (MVA2), TRIP13 (MVA3), CENATAC (MVA4), MAD1L1 (MVA7), and BUB1
  and BUB3 (Orphanet MVA group), and builds constraint and NMD maps for each
  (`findings/genes/`). CENATAC and MAD1L1 had been missed by the hand-made list.

## 9. Status and possible extensions

- Done: report (`docs/smarin00_track2_report.md` and `.pdf`), methods form,
  pitch video.
- Possible extension: physics- or ML-based ΔΔG (e.g. ThermoMPNN) as a further
  line of evidence, calibrated against the same six alleles.

## References

- Hanks S. et al. Constitutional aneuploidy and cancer predisposition caused
  by biallelic mutations in BUB1B. *Nat Genet* 36, 1159–1161 (2004).
- Suijkerbuijk S.J.E. et al. Molecular causes for BUBR1 dysfunction in the
  human cancer predisposition syndrome mosaic variegated aneuploidy.
  *Cancer Res* 70, 4891–4900 (2010).
- Suijkerbuijk S.J.E. et al. The vertebrate mitotic checkpoint protein BUBR1
  is an unusual pseudokinase. *Dev Cell* 22, 1321–1329 (2012).
- Choi E. et al. *EMBO J* 28, 2077–2089 (2009). doi:10.1038/emboj.2009.123
- Park I. et al. *J Cell Biol* 202, 295–309 (2013). doi:10.1083/jcb.201210099
- Park I. et al. *FEBS J* 284, 4035–4050 (2017). doi:10.1111/febs.14286
- BubR1 TPR domain supports the mitotic checkpoint by promoting MCC formation
  and MCC-APC/C interaction. *Commun Biol* (2026),
  doi:10.1038/s42003-026-10384-9.
