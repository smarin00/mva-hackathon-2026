# Track 2: Drug candidates

Status: 19 Oct 2026. Builds on `track2_conclusions.md` (mechanism) and
`findings/track2/opentargets/ot_*.tsv` (automated layer, Open Targets Platform API 26.9.0,
data 26.09).

## Therapeutic premise

The child's BubR1 comes entirely from the N1002K allele (L737\* is
NMD-sensitive). N1002K most likely destabilises the protein, as the
neighbouring MVA allele L1012P does. MVA1 is a disorder of BubR1 *quantity*
with a steep threshold (~6% of normal: widespread missegregation; ~13%:
virtually none). **The goal is to raise the steady-state level of an existing,
functional-but-unstable protein by roughly two-fold.** There is no enzymatic
activity to restore (the kinase domain is a pseudokinase).

Two routes to more BubR1 protein, and one to its downstream consequences:

| Route | Lever | Evidence for BubR1 specifically |
|---|---|---|
| A. Reduce BubR1 degradation | NAD⁺ → SIRT2 → deacetylation of BubR1 K668 | Direct, in vivo (North 2014) |
| B. Improve folding of the mutant | Heat-shock response / HSP70 chaperones | Mechanistic extrapolation |
| C. Limit downstream damage | Clear senescent cells | Genetic, in BubR1 progeroid mice (Baker 2011) |

## Ranked candidates

The brief asks for **existing, approved medications**. Each candidate below is
an approved drug; non-approved alternatives with the same mechanism are noted
where they have a better safety or pediatric record.

### 1. Niacin (nicotinic acid), NAD⁺ repletion, lead candidate

- **Approved drug.** FDA-approved prescription medicine (dyslipidaemia).
- **NAD⁺ effect in humans.** 750–1,000 mg/day raised blood NAD⁺ up to 8-fold
  and restored muscle NAD⁺ in adult mitochondrial myopathy (Pirinen et al.,
  Cell Metab 2020).
- **Why niacin.** It enters NAD⁺ synthesis through NAPRT (Preiss-Handler
  pathway), bypassing NAMPT, and unlike nicotinamide it does not feed back to
  inhibit sirtuins.
- **Risks.** Flushing, hepatotoxicity at high dose (liver monitoring), an
  "anaemia tendency" in Pirinen 2020; no pediatric NAD⁺-repletion data.
- **Non-approved alternative:** nicotinamide riboside (dietary supplement),
  better tolerated, with pediatric precedent (below).

The mechanism and evidence for the NAD⁺ route:

- **Mechanism.** BubR1 levels fall when NAD⁺ and SIRT2 activity fall. SIRT2
  deacetylates BubR1 at K668, which CBP acetylates, and increased SIRT2 raises
  BubR1 abundance in vivo. The NAD⁺ precursor NMN raised BubR1 levels in mice.
  Raising SIRT2 **extended the median lifespan of BubR1 hypomorphic mice**,
  the closest animal model of MVA's dosage loss (North et al., EMBO J 2014).
- **Fit to this child.** N1002K is full-length and retains K668, so the
  pathway still has its substrate.
- **Pediatric precedent.** NR 25 mg/kg/day for 4 months in 24 patients with
  ataxia-telangiectasia, 17 of them under 18 (Movement Disorders 2021): a
  pediatric genome-instability disorder.
- **NR status.** Dietary supplement (GRAS in the US), not an approved
  medicine. Open Targets reports "APPROVAL" for NR; this does not reflect
  regulatory approval.
- **Why NR, not nicotinamide.** Nicotinamide is the cheaper NAD⁺ precursor,
  but it is also the product-feedback inhibitor of sirtuins and inhibits them
  at high concentrations. NR and NMN enter NAD⁺ synthesis downstream of that
  problem.
- **Gap.** SIRT2-mediated stabilisation was shown for wild-type BubR1; it is
  untested whether it rescues a misfolding-prone mutant.
- **Second acetylation site, opposite sign.** BubR1 also carries K250, whose
  acetylation by PCAF *protects* it from APC/C-mediated degradation in
  mitosis (Choi 2009). SIRT2 deacetylated K250 in vitro (Suematsu 2014, as
  cited in Park 2017), which in principle works against BubR1. Park 2017 found
  no mitosis-specific SIRT2–BubR1 binding and suggests SIRT2 acts mainly in
  interphase; the in vivo net effect of more SIRT2 was more BubR1 and longer
  life (North 2014). The validation plan measures BubR1 through mitosis, not
  only at steady state, for this reason.
- **Measurable.** BubR1 protein in patient lymphocytes; whole-blood NAD⁺.

### 2. Arimoclomol (heat-shock response amplifier)

- **Mechanism.** Amplifies HSF1-driven induction of HSP70-family chaperones
  under cellular stress, supporting folding of misfolding-prone proteins.
  N1002K is predicted to destabilise BubR1, and its neighbour L1012P is a
  verified chaperone-dependent client.
- **Status.** FDA-approved September 2024 (Miplyffa) for Niemann-Pick disease
  type C, **adults and children aged 2 and older**, oral. **Not approved in
  the EU:** the EMA's CHMP adopted a negative opinion in July 2026, and a
  re-examination is pending. Open Targets still lists it as phase 3.
- **Gap.** No BubR1-specific data; this is a mechanism-class extrapolation.
  The HSF1 response also induces HSP90, which should help an HSP90 client,
  but that is a prediction, not an observation.

### 3. Senolytic: dasatinib + quercetin (tier 3, downstream)

- **Mechanism.** Clearing p16Ink4a-positive senescent cells in BubR1 progeroid
  mice delayed ageing phenotypes in eye, muscle and fat, even when started
  late (Baker et al., Nature 2011).
- **Approved drug.** Dasatinib is FDA-approved, including for children with
  Ph+ ALL. Quercetin is a dietary flavonoid. Intermittent D+Q dosing is the
  senolytic regimen studied in humans.
- **Risks.** Dasatinib is a multi-kinase oncology drug (myelosuppression,
  pleural effusion). Fisetin is a better-tolerated, non-approved alternative;
  navitoclax is excluded (dose-limiting thrombocytopenia).
- **Gap.** Baker 2011 used *genetic* clearance, not drugs. This treats
  consequences, not cause. Long-term and lowest priority.

## Considered, not ranked: HDAC inhibition (valproic acid)

HDAC2 and HDAC3 deacetylate BubR1 at K250, so a class I HDAC inhibitor would
keep BubR1 acetylated. Valproic acid is an approved class I HDAC inhibitor
with long pediatric use (epilepsy), so it is the obvious approved drug for
this route. We read the three primary papers in full; they argue against
ranking it.

- **It acts on timing in mitosis, not on dose.** K250 is acetylated only in
  prometaphase, and HDAC2/3 bind BubR1 only in mitosis, not in interphase
  (Choi 2009; Park 2017). Acetylation decides *when* BubR1 is destroyed at
  the end of mitosis. An acetylation-deficient allele left overall BubR1
  levels unchanged in mouse tissues (Park 2013). This child's deficit is
  BubR1 *quantity*, set largely by turnover of an unstable protein; HDAC
  inhibition does not address that.
- **The direct effect is mitotic arrest.** Keeping BubR1 acetylated (the
  K250Q mimic) held cells in metaphase about 45 min longer than wild type
  (Park 2017); the pan-HDAC inhibitor TSA delayed mitotic exit by at least
  5 h (Choi 2009). HDAC inhibitors act against cancer cells partly through
  mitotic arrest (Park 2017). In MVA, where mitosis is too short, a longer
  checkpoint could in principle help; but whether it helps cells with too
  little BubR1, or simply arrests and kills dividing cells (marrow, gut), is
  untested.
- **Where it could help.** Acetylated BubR1 is needed to recruit PP2A-B56 for
  stable kinetochore attachment, and the acetylation-deficient allele caused
  tumours and missegregation in mice (Park 2013). Raising the acetylated
  fraction of the little BubR1 this child has is a coherent hypothesis.
- **Pharmacology.** Valproate is a weak, broad class I inhibitor acting at
  millimolar concentrations; its effect on BubR1 K250 at therapeutic levels is
  unknown. Known risks include hepatotoxicity (highest under age 2) and
  pancreatitis.

Verdict: a mechanistically coherent but double-edged hypothesis. It is
included in the validation plan as an exploratory arm (mitotic timing and
missegregation in patient cells), not as a ranked candidate.

## Rejected at curation

| Agent | Why it appeared | Why rejected |
|---|---|---|
| PRI-724 | CREBBP inhibitor (direction-consistent) | Blocks the CBP–β-catenin interaction, not CBP acetyltransferase activity, so it cannot affect BubR1 K668 acetylation. Also a Wnt modulator. |
| Navitoclax, obatoclax | BCL-xL inhibitors (senolytic direction) | Oncology agents; thrombocytopenia; not justified in a child when fisetin exists. |
| Nicotinamide (high dose) | NAD⁺ precursor | Sirtuin inhibitor at high concentrations (see above). |
| Proteasome inhibitors | Rescued mutant BubR1 levels in vitro | Systemic toxicity; mechanism proof, not a treatment. |

## Do-no-harm list

Drugs that act *against* the therapeutic direction. Several are in oncology
pipelines, which matters for a child with MVA-associated cancer
predisposition (rhabdomyosarcoma, Wilms tumour, leukaemia).

| Class | Count in Open Targets | Mechanism of harm | Examples |
|---|---|---|---|
| HSP90 inhibitors | 10 | Mutant BubR1 is an HSP90 client; inhibition depletes it (Suijkerbuijk 2010) | ganetespib, tanespimycin, retaspimycin, luminespib, onalespib |
| NAMPT inhibitors | 3 | Deplete NAD⁺, opposing candidate 1 | daporinad (trialled in childhood leukaemia), teglarinad |
| Aurora B inhibitors | 18 | Impair error correction in an already weakened checkpoint | barasertib, ibcasertib, tozasertib, danusertib |
| MPS1 / TTK inhibitors | 2 | Weaken the spindle assembly checkpoint directly | BAY-1161909, BAY-1217389 |
| SIRT2 inhibitors | 0 in clinic | Would block the K668 deacetylation route | AGK2 (research tool) |

## Proposed validation (before any clinical use)

In patient-derived cells (lymphoblastoid line or fibroblasts):

1. Baseline BubR1 protein level vs. healthy controls (western blot); expected
   well below 50%.
2. Treat with nicotinic acid (and NR for comparison), and separately arimoclomol: measure BubR1 level,
   nocodazole-induced checkpoint arrest, and missegregation rate (live imaging
   or FISH aneuploidy count).
3. HSP90 inhibitor as a positive control for further depletion, which also
   tests the do-no-harm prediction directly.
4. Cycloheximide chase to show any rise in level comes from slower turnover.
5. Exploratory: valproic acid at therapeutic concentrations (0.3–0.7 mM),
   reading out acetyl-K250 BubR1, time from nuclear-envelope breakdown to
   anaphase, and missegregation. A longer mitosis without fewer errors, or
   mitotic death, would rule it out.

Success criterion: BubR1 raised from below to above ~13% of normal, with a
measurable fall in missegregation.

## Automated layer: how candidates were screened

`bin/opentargets_candidates.py` queries 13 pathway genes (`config/bub1b.toml`), each with a
literature-justified desired direction, against Open Targets known drugs and
clinical candidates, and classifies every drug by its ChEMBL action type.

- Result: no drug acts directly on BUB1B, SIRT2 or HSF1. Leads 1 and 2 enter
  through named-drug lookup and literature curation, not target matching.
- Database limitations found: stale regulatory status (arimoclomol), and
  "APPROVAL" flags on supplements (NR, quercetin, fisetin). Regulatory status
  in this report is taken from primary sources.

## References

- North B.J. et al. SIRT2 induces the checkpoint kinase BubR1 to increase
  lifespan. *EMBO J* 33, 1438–1453 (2014). doi:10.15252/embj.201386907
- Baker D.J. et al. Clearance of p16Ink4a-positive senescent cells delays
  ageing-associated disorders. *Nature* 479, 232–236 (2011).
  doi:10.1038/nature10600
- Suijkerbuijk S.J.E. et al. *Cancer Res* 70, 4891–4900 (2010).
- Choi E. et al. BubR1 acetylation at prometaphase is required for modulating
  APC/C activity and timing of mitosis. *EMBO J* 28, 2077–2089 (2009).
  doi:10.1038/emboj.2009.123
- Park I. et al. Loss of BubR1 acetylation causes defects in spindle assembly
  checkpoint signaling and promotes tumor formation. *J Cell Biol* 202,
  295–309 (2013). doi:10.1083/jcb.201210099
- Park I. et al. HDAC2/3 binding and deacetylation of BubR1 initiates spindle
  assembly checkpoint silencing. *FEBS J* 284, 4035–4050 (2017).
  doi:10.1111/febs.14286
- Suematsu T. et al. Deacetylation of the mitotic checkpoint protein BubR1 at
  lysine 250 by SIRT2 and subsequent effects on BubR1 degradation during the
  prometaphase/anaphase transition. *Biochem Biophys Res Commun* (2014).
  Cited via Park 2017; not read in full.
- Pirinen E. et al. Niacin cures systemic NAD+ deficiency and improves muscle
  performance in adult-onset mitochondrial myopathy. *Cell Metab* 31,
  1078–1090 (2020). doi:10.1016/j.cmet.2020.04.008
- Zevra Therapeutics. Negative CHMP opinion for arimoclomol (24 Jul 2026).
- Nicotinamide riboside in ataxia-telangiectasia. *Movement Disorders* (2021).
  doi:10.1002/mds.28788
- FDA. FDA approves first treatment for Niemann-Pick disease type C
  (20 Sep 2024).
