# Placements

**Related:** [`../PRODUCT.md`](../PRODUCT.md) — the principles this design applies;
[`../../GLOSSARY.md`](../../GLOSSARY.md) — assembly, alignment, alternate locus, patch.

## Overview

An assembly is not one sequence per chromosome. GRCh38.p14 has 24 chromosomes, a mitochondrion, 355 `NT_` scaffolds and
325 `NW_` patches, and NCBI aligns every RefSeq transcript to every one of them that carries its gene. The store keeps
every one of those placements, ordered so that a caller who names no sequence gets the chromosome. The CDS is taken from
the placements that state both its ends; where they disagree, the record's own protein arbitrates. A transcript none of
whose placements states its CDS is marked so, and the provider refuses to model it rather than serve it as non-coding. A
gene's features across sequences are one bundle, joined by GeneID.

The reason is the third principle in [`PRODUCT.md`](../PRODUCT.md): an answer or an error, never a guess. A gene the
primary assembly's haplotype lacks has no chromosome position, and the honest answer to "where is this variant on
GRCh38" is its position on the alternate locus that carries the gene. A transcript whose 5' end a fix patch supplies has
two honest answers, and which one a caller wants depends on what they are doing with it.

## Background

The GRC releases three kinds of sequence beside the chromosomes. **Alternate loci** (`NT_`, in `ALT_REF_LOCI_*`) are
other haplotypes of regions too variable for one sequence: the MHC, the KIR cluster. **Unlocalized and unplaced
scaffolds** are also `NT_`: sequence known to belong to a chromosome, or to the genome, but not where. **Patches**
(`NW_`) come between assembly versions: a fix patch corrects the primary assembly and is folded into the next major
release; a novel patch adds alternate sequence and becomes an alternate locus.

NCBI annotates all of them. A gene present on a chromosome and on a patch is two gene features with a suffixed feature
id and one GeneID; each of its transcripts is placed once per sequence, and each placement states the CDS in that
sequence's coordinates. Where the record's coding sequence runs past what the sequence carries, the CDS rows are marked
partial and the open end is stated — RYBP's `NM_012234.7` on chromosome 3, whose first 213 bases the primary assembly
lacks. NCBI marks every row partial when *any* end of *any* row is open, which includes the two rows either side of an
internal break such as a frameshift, so the flag alone does not say whether the CDS's own bounds are stated.

Some genes are on no chromosome at all: GSTT1, where the reference chromosome 22 carries the deletion allele; HLA-DRB3
and HLA-DRB4, which the reference MHC haplotype lacks; the KIR genes absent from the reference KIR haplotype. On
GRCh38.p14, 8,029 transcripts are placed only on scaffolds or patches. Nineteen MANE Select transcripts are on alternate
loci and 45 on patches.

## Non-goals

- **No projection between placements.** A variant on an alternate locus and one on the chromosome are on different
  sequences; relating them is liftover, which [`PRODUCT.md`](../PRODUCT.md) puts outside this package.
- **No choice of placement on the caller's behalf beyond the ordering.** The provider serves the placement named, or the
  first; which sequence a study means is the study's decision.
- **No CDS the annotation does not state.** The CDS could be recovered from the protein sequence, or read from the
  record's GenBank entry; the builder does neither. See Alternatives.

## Design

### Every placement is bundled

A transcript's `alignments` hold one entry per sequence NCBI aligns it to, chromosomes first, then scaffolds, then
patches, each by accession. The order is a promise of the at-rest format
([`bundle.proto`](../../proto/weaver_data_provider/v1/bundle.proto)): the provider's `get_transcript` with no sequence
named returns the first, so a transcript on a chromosome is modelled there by default, and one on scaffolds only is
modelled on the lowest-numbered of them. A caller who wants another names it, and a sequence the transcript is not
placed on is an error.

Consequences: a store's interval tables have one table per sequence, so a position query on `NT_187633.1` finds GSTT1;
the genome must be cut from the assembly FASTA that includes the scaffolds and patches; a `c.` variant on an
alternate-locus-only gene projects to `g.` on that locus and to nothing on a chromosome.

### The CDS comes from the placements that state it whole

Each placement's CDS rows are read with the openness of each end. A placement determines the CDS when it has CDS rows
and neither the lowest row's start nor the highest row's end is open. Every determining placement with an alignment
projects its bounds through that alignment to transcript indices. They usually agree. Where they do not — NIPA2's
`NM_001008860.3` on GRCh38.p14 projects its CDS end one base earlier through a patch whose alignment has indels in the
last exon than through chromosome 15 — an alignment or an annotation is wrong, and the record says which: the CDS is the
projection over which the record, read as whole codons in the standard genetic code, encodes its own protein and then a
stop. The other placements keep their alignments, lose their say over the CDS, and are counted in the build's report.
The reading is deliberately strict: a selenocysteine, a non-AUG initiator or an ambiguous base reads differently, so a
record with one of those and disagreeing placements fails the build, naming translation as the reason, rather than being
arbitrated by a looser rule. The build also fails when the record names no protein, the protein set lacks it, or the
record encodes it over none of the projections, since a bundle would then be guessing.

A coding transcript with no placement that both has an alignment and states the CDS whole — 32 of 136,269 on GRCh38.p14,
all RefSeq models whose 5' end no sequence carries — is bundled with its sequence, protein and alignments, no CDS, and
`cds_undetermined` set. The flag is a fact about the bundle, not the annotation: a placement the annotation states a
whole CDS on but NCBI published no alignment for projects nothing. The provider refuses to model such a transcript,
naming the reason, since serving it with no CDS would make weaver read every `c.` name on it as non-coding: a wrong
answer rather than a missing one. The refusal is a change an older reader could not know to make, so it comes with a new
format version, which the genome catalogue shares: a genome in a bucket is rebuilt with its store although its layout
did not change.

### A gene is one bundle across sequences

Gene features are grouped by GeneID, which every NCBI gene feature carries and which the copies of one gene on a
chromosome, a patch and an alternate locus share. A feature without one fails the build.

## Alternatives considered

- **Chromosomes only.** Simple, and wrong for 8,029 transcripts: an alternate-locus-only gene is not in the store, and a
  variant in a fix-patch gene has only the soft-clipped chromosome placement to be projected through.
- **First placement's CDS wins.** The first might be the partial one, and the check that the placements agree is lost.
- **Fail the build when placements disagree.** One transcript on GRCh38.p14 disagrees, so this refuses every release for
  one patch alignment's indel, when the record itself can settle it.
- **Prefer the chromosome placement when placements disagree.** Deterministic and usually right, but a preference is a
  guess where the record's protein is a fact; the two coincide on GRCh38.p14, and the rule that reads the record stays
  right when they do not.
- **Arbitrate with a tolerant reading** — `U` for an internal stop, `M` for any initiator, a wildcard for an ambiguous
  base. Covers RefSeq's real exceptions, and accepts a projection that is not the record's. A strict reading and a build
  failure name the case instead.
- **Fail the build on a coding transcript with no determined CDS.** The principled default for a missing input, and what
  a coding transcript with no alignment gets. Rejected because 32 records would refuse the whole release, and the state
  is a fact about the publisher's data, not a broken input; making it representable and refusing at the provider keeps
  both the build and the honesty.
- **Serve such a transcript as non-coding.** Rejected: a wrong answer, indistinguishable at rest from a non-coding
  transcript.
- **Read the CDS from the record** (NCBI's `rna.gbff`), which states it in transcript coordinates and would remove the
  undetermined class and the per-placement projection. Not taken: a further multi-gigabyte input for 32 records, and the
  projection through the alignment is also a check that the alignment and the annotation agree.
- **Derive the CDS from the protein**, as the open reading frame whose translation is the protein. A derivation, with
  its own edge cases (selenoproteins, non-AUG initiators, ribosomal slippage), which the first principle rules out; the
  arbitration above reads the protein only to choose between bounds the annotation stated.
- **Classify `NT_` scaffolds** into alternate loci and unlocalized scaffolds from the assembly report, for the report
  and the ordering. Not done: NCBI's annotation does not distinguish them, both sort after the chromosomes, and the
  order among them is by accession either way.

## Open questions

- Whether a caller of `get_transcript` with no sequence named should be told, in the model returned, that other
  placements exist. weaver's `TranscriptData` has no field for it today.
- How Ensembl placements join this design: Ensembl transcripts are genome-derived, so their placements are their exon
  coordinates, and a MANE Select transcript on a fix patch has the RefSeq partner's NCBI alignment to the chromosome
  available for copying. To be decided with the Ensembl builder.
- A region query lists a transcript with an undetermined CDS among its answers, and a caller that then models each one
  meets the refusal for it. Whether the listing should say so first is open.
