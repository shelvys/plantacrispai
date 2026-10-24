# Scientific Methods and Attribution

The published methods PlantaCRISPAI implements or adapts, and the limits of what
it claims. Implementation detail is in [ARCHITECTURE.md](ARCHITECTURE.md).

## 1 · Sequence alignment

Local alignment uses the **Smith–Waterman** algorithm with affine-like gap
handling, anchored to the user-supplied reference amplicon.

> Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular
> subsequences. *Journal of Molecular Biology*, 147(1), 195–197.
> https://doi.org/10.1016/0022-2836(81)90087-5

Reads are auto-oriented before alignment. Reverse-primer reads are detected by
comparing forward and reverse-complement alignment scores and corrected, since an
unoriented reverse read aligns poorly and is readily miscalled as heavily edited.

## 2 · Decomposition-based editing quantification

The editing score is **ICE-equivalent**: it reproduces the interpretation of
decomposition-based Sanger analysis — indel size, post-cut trace degradation and
inferred allele composition — from alignment-derived quantities rather than by
deconvolving the fluorescence trace itself.

> Conant, D., Hsiau, T., Rossi, N., Oki, J., Maures, T., Waite, K., Yang, J.,
> Joshi, S., Kelso, R., Holden, K., Enzmann, B. L., & Stoner, R. (2022).
> Inference of CRISPR edits from Sanger trace data. *The CRISPR Journal*, 5(1),
> 123–130. https://doi.org/10.1089/crispr.2021.0113

> Brinkman, E. K., Chen, T., Amendola, M., & van Steensel, B. (2014). Easy
> quantitative assessment of genome editing by sequence trace decomposition.
> *Nucleic Acids Research*, 42(22), e168. https://doi.org/10.1093/nar/gku936

> Pinello, L., Canver, M. C., Hoban, M. D., Orkin, S. H., Kohn, D. B., Bauer,
> D. E., & Yuan, G.-C. (2016). Analyzing CRISPR genome-editing experiments with
> CRISPResso. *Nature Biotechnology*, 34(7), 695–697.
> https://doi.org/10.1038/nbt.3583

> [!IMPORTANT]
> PlantaCRISPAI does not reimplement ICE, TIDE or CRISPResso, and its scores are
> not numerically interchangeable with theirs. The ICE-equivalent score is
> calibrated to agree in *classification*, not in value. Where an exact ICE or
> TIDE number is required — for a manuscript figure, say — run those tools.

## 3 · Guide RNA efficiency

The efficiency model uses position-specific nucleotide weights derived from
large-scale screens, with plant-specific corrections.

> Doench, J. G., Fusi, N., Sullender, M., Hegde, M., Vaimberg, E. W., Donovan,
> K. F., Smith, I., Tothova, Z., Wilen, C., Orchard, R., Virgin, H. W.,
> Listgarten, J., & Root, D. E. (2016). Optimized sgRNA design to maximize
> activity and minimize off-target effects of CRISPR-Cas9. *Nature
> Biotechnology*, 34(2), 184–191. https://doi.org/10.1038/nbt.3437

> Wang, T., Wei, J. J., Sabatini, D. M., & Lander, E. S. (2014). Genetic screens
> in human cells using the CRISPR-Cas9 system. *Science*, 343(6166), 80–84.
> https://doi.org/10.1126/science.1246981

Plant-specific feature weighting draws on:

> Xu, H., Xiao, T., Chen, C.-H., Li, W., Meyer, C. A., Wu, Q., Wu, D., Cong, L.,
> Zhang, F., Liu, J. S., Brown, M., & Liu, X. S. (2015). Sequence determinants of
> improved CRISPR sgRNA design. *Genome Research*, 25(8), 1147–1157.
> https://doi.org/10.1101/gr.191452.115

> Endo, M., Mikami, M., & Toki, S. (2015). Multigene knockout utilizing
> off-target mutations of the CRISPR/Cas9 system in rice. *Plant and Cell
> Physiology*, 56(1), 41–47. https://doi.org/10.1093/pcp/pcu154

> Tang, X., Liu, G., Zhou, J., Ren, Q., You, Q., Tian, L., Xin, X., Zhong, Z.,
> Liu, B., Zheng, X., Zhang, D., Malzahn, A., Gong, Z., Qi, Y., Zhang, T., &
> Zhang, Y. (2018). A large-scale whole-genome sequencing analysis reveals highly
> specific genome editing by both Cas9 and Cpf1 nucleases in rice. *Genome
> Biology*, 19, 84. https://doi.org/10.1186/s13059-018-1458-5

> [!NOTE]
> The Doench position weights were derived in human cells. Their transfer to
> plant systems is an approximation. Scores are useful for *ranking* candidate
> guides within one target region; they are not calibrated absolute efficiency
> predictions for any plant species.

## 4 · Thermodynamics

Nearest-neighbour model. RNA:DNA parameters for guide:target duplex stability;
DNA:DNA parameters for primer melting temperature.

> Sugimoto, N., Nakano, S., Katoh, M., Matsumura, A., Nakamuta, H., Ohmichi, T.,
> Yoneyama, M., & Sasaki, M. (1995). Thermodynamic parameters to predict
> stability of RNA/DNA hybrid duplexes. *Biochemistry*, 34(35), 11211–11216.
> https://doi.org/10.1021/bi00035a029

> SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide
> DNA nearest-neighbor thermodynamics. *PNAS*, 95(4), 1460–1465.
> https://doi.org/10.1073/pnas.95.4.1460

> Freier, S. M., Kierzek, R., Jaeger, J. A., Sugimoto, N., Caruthers, M. H.,
> Neilson, T., & Turner, D. H. (1986). Improved free-energy parameters for
> predictions of RNA duplex stability. *PNAS*, 83(24), 9373–9377.
> https://doi.org/10.1073/pnas.83.24.9373

Duplex ΔG, seed-region ΔG and self-folding penalty are computed at 22 °C, a
representative plant transformation temperature.

## 5 · CAPS and dCAPS markers

Cleaved Amplified Polymorphic Sequence markers detect an edit by whether a
restriction site survives it. Derived CAPS extends this by introducing a
deliberate primer mismatch that completes a restriction site in one allele only.

> Konieczny, A., & Ausubel, F. M. (1993). A procedure for mapping *Arabidopsis*
> mutations using co-dominant ecotype-specific PCR-based markers. *The Plant
> Journal*, 4(2), 403–410. https://doi.org/10.1046/j.1365-313X.1993.04020403.x

> Neff, M. M., Neff, J. D., Chory, J., & Pepper, A. E. (1998). dCAPS, a simple
> technique for the genetic analysis of single nucleotide polymorphisms.
> *The Plant Journal*, 14(3), 387–392.
> https://doi.org/10.1046/j.1365-313X.1998.00124.x

Enzyme recognition sequences, cut positions, overhang types and incubation
temperatures follow supplier catalogue data (New England Biolabs; Thermo Fisher
Scientific). **Always confirm against the current supplier datasheet before
ordering or digesting** — catalogue entries change.

## 6 · File format

`.ab1` files follow the **ABIF** (Applied Biosystems Instrument File) container
specification. The parser reads the format for interoperability only.

> Applied Biosystems. *ABIF File Format Specification and Sample Data File*.
> Thermo Fisher Scientific.

## 7 · Machine learning

A Random Forest classifier over 45 alignment- and quality-derived features is
ensembled with rule-based and ICE-equivalent scores. The guide designer uses a
separate 28-feature model.

> Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5–32.
> https://doi.org/10.1023/A:1010933404324

### Validation status — read this before citing a number

The classifier was developed and evaluated on an *OsAAP3* rice T1 panel
independently genotyped by ICE and TIDE, spanning wild-type, mosaic,
heterozygous and homozygous classes.

> [!WARNING]
> **The reported accuracy is inconsistent across the current codebase.** Three
> different figures appear:
>
> | Location | Figure |
> |---|---|
> | `pages/LandingPage.js` | 89.5 % on 23 rice T1 samples |
> | `engine/analysisEngine.js` | 18/19 correct (94.7 %) |
> | `engine/demoData.js` | 17/20 confirmed samples correct |
>
> These appear to reference different snapshots of the validation set taken
> during development. **They must be reconciled against the source validation
> table before any manuscript, thesis or regulatory submission.** Until then,
> treat the figure as provisional and cite the underlying table, not the app.

Further limitations of the validation as it stands:

- **Single gene, single species, single generation.** All validation samples are
  *OsAAP3* rice T1. Performance on other genes, species, generations or ploidy
  levels is untested.
- **No held-out test set.** Thresholds were tuned on the same panel used to
  report accuracy, so the figure is optimistic. An independent test set is
  required for a defensible performance claim.
- **Small n.** With roughly 20 samples, a single reclassification moves accuracy
  by about 5 %. Confidence intervals are wide.
- **Ground truth is itself inferential.** ICE and TIDE are computational calls,
  not sequencing of cloned alleles. Where they disagree — as recorded for sample
  223 — there is no gold standard.

The honest statement of current status: PlantaCRISPAI agrees closely with ICE and
TIDE on a small rice T1 panel. Broader validation is ongoing.

## 8 · Off-target assessment

Off-target risk is scored **heuristically** from seed-region composition and
sequence features. There is no genome-wide alignment, because a genome index
cannot be shipped to a browser.

> [!CAUTION]
> Do not treat the off-target score as a genome-wide specificity assessment.
> Check candidate guides against the target genome with a dedicated tool before
> ordering oligos.

> Concordet, J.-P., & Haeussler, M. (2018). CRISPOR: intuitive guide selection
> for CRISPR/Cas9 genome editing experiments and screens. *Nucleic Acids
> Research*, 46(W1), W242–W245. https://doi.org/10.1093/nar/gky354

> Bae, S., Park, J., & Kim, J.-S. (2014). Cas-OFFinder: a fast and versatile
> algorithm that searches for potential off-target sites of Cas9 RNA-guided
> endonucleases. *Bioinformatics*, 30(10), 1473–1475.
> https://doi.org/10.1093/bioinformatics/btu048

## 9 · Regulatory context

Generation-aware guidance reflects India's 2022 regulatory framework exempting
SDN-1 and SDN-2 genome-edited plants free of foreign DNA from GMO rules.

> Ministry of Environment, Forest and Climate Change, Government of India (2022).
> *Guidelines for the Safety Assessment of Genome Edited Plants*.

> Department of Biotechnology, Government of India (2022). *Office Memorandum
> on the exemption of SDN-1 and SDN-2 genome edited plants.*

> [!NOTE]
> Regulatory guidance changes. Confirm current requirements with your
> institutional biosafety committee. Nothing in this software is regulatory
> advice.

## Disclaimer

Citations credit the authors of methods this software implements or adapts.
They do not imply that those authors endorse, reviewed or validated
PlantaCRISPAI.
