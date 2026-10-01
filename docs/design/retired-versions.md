# Retired transcript versions

**Related:** [`../PRODUCT.md`](../PRODUCT.md) — the principles this design applies; [`placements.md`](placements.md) —
how a transcript's placements and CDS are bundled, which a retired version shares;
[`../../GLOSSARY.md`](../../GLOSSARY.md) — annotation release, shard, store.

## Overview

A variant name cites a transcript version, and the version it cites outlives the version's currency: ClinVar's
expressions name `NM_000059.3` for years after NCBI has moved to `.4`. An annotation release names only the current
versions, so a store built from releases alone has no model for the names most of the world's variant records carry.
NCBI publishes a remedy beside each GRCh38 annotation release: a historical set of every replaced and suppressed `NM_`
and `NR_` version with the alignment it last had, as the same kinds of file as a release. The store takes that set as a
shard of its own, stacked under the current release's shard, so a retired version is served when named and the current
release's model wins for a version both hold. Each version carries NCBI's current status for it — current, superseded or
suppressed — read not from its record, which does not know when it was retired, but from a table fetched from Entrez and
given to the build. A retired version carries no tags and no MANE partner, and a region's answer lists current versions
only: what the publisher recommends today is a property of today's release, not of a model that was once current.

## Background

A RefSeq transcript version is immutable: a change to the sequence is a new version, and the old one is *replaced*; a
record withdrawn altogether is *suppressed*. The old version remains in GenBank and in Entrez, where its summary states
its status and, for a replaced one, its successor — which may be an unversioned accession, or another accession the
record was merged into. The record's own COMMENT states a replacement only when the replacement happened before the
record was last written, and never states a suppression.

NCBI's historical set (`RefSeq_historical_alignments` under each GRCh38 annotation release) is cumulative: every known
version ever aligned to GRCh38, each with the alignment from the most recent release that carried it, supplemented with
older versions never annotated there, and filtered of alignments far from the gene's expected location. It is anchored
on one release — RS_2023_03 at the time of writing — and about half its 189,000 versions are retired. It comes as an
annotation GFF3 generated from the alignments, the alignment BAM with the records' sequences, and the records as a
GenBank flat file, but no sequence FASTA sets: the records carry each transcript's sequence and, on its CDS feature, the
protein id and translation. There is no such set for GRCh37.

The provider resolves a named version from the newest shard that holds it and lists a region's transcripts from every
shard ([`placements.md`](placements.md) has the store's shape).

## Non-goals

- **No resolution of a retired version to its successor.** Entrez names the successor, and a caller wanting the current
  model's answer for an old name can ask for it; which version a study means is the study's decision.
- **No GRCh37 historical shard.** NCBI publishes none; building one from older release files would be a derived set.
- **No status for an Ensembl transcript.** Ensembl archives releases rather than retiring versions; a retired Ensembl
  model is in an older release's shard, if a store stacks one.

## Design

### The historical set is a shard like any release's

The RefSeq builder reads the set's GFF3, BAM and GenBank files as it reads a release's, with two differences the
builder's `historical` command expresses. The transcript and protein sequences come from the GenBank records, since the
set has no FASTA sets: the ORIGIN section is the record, the CDS feature's translation is the protein, and a record read
for its sequences without either fails the build. And a status table is an input, with every version named in the set,
so that each transcript's status is a fact stated by the publisher rather than one the builder inferred.

The shard is stacked under the current release's shard. The store's rule that the newest shard holding a version answers
for it then gives a version in both — the roughly half of the set that was current at its anchor — the current release's
model, and a version only the historical set holds its last alignment. Nothing in the store or the provider changes for
this: the shard's immutability and the stacking order already carry it.

### Status comes from Entrez, through a table

The status of a version is a fact about the publisher's catalogue today. The record cannot state it, as Background says,
and the current annotation cannot either: a version absent from it may be live and merely unannotated (`NM_201563.5` on
RS_2024_08), and a version whose accession is absent may have been replaced rather than suppressed (`NM_000028.2`,
replaced by a `.3` that post-dates the release). Entrez states it, per version: live, replaced with a successor, or
suppressed.

The build reads local files only, so the Entrez answers are fetched beforehand by a script in the repository into a
table of one row per version, and the table is an input to the build, recorded with the shard's other inputs by digest.
A version with no row, or a row with a status word that is not one of Entrez's, fails the build. The table is a
snapshot: a version live when the table was fetched and replaced since is current in the shard until the shard is cut
again. The alternative, a build that queries Entrez itself, would make a shard's contents depend on the network at the
moment of the build and leave no input to record.

### A retired version carries no recommendation

The annotation's `tag` attribute and the MANE summary name what the publisher recommends: RefSeq Select, MANE Select. A
tag on a retired version in the set's GFF3 says what was recommended while that version was current, and the MANE
summary names current versions only. A retired version is bundled with neither, so a caller reading tags reads only
recommendations that hold.

For the same reason the provider's region lookup lists current transcripts only. A retired version is served when named
— that is the shard's purpose — but a caller asking "which transcripts are here" is asking about today's annotation, and
the historical shard would otherwise answer with a decade of replaced models for every gene.

## Alternatives considered

- **Infer status from the current annotation**: a version it names is current, one whose accession it names under a
  newer version is superseded, the rest suppressed. Wrong in both directions on real data, as Design says.

- **Take the record's COMMENT as the status.** It names 97,821 of the set's versions as replaced and is silent on the
  other 4,356 retired ones, so a shard built from it would hold suppressed records as current.

- **Query Entrez from the build.** Honest about the status, but a build whose inputs are not all files cannot be
  recorded or repeated.

- **A shard per past annotation release**, stacked oldest first, instead of NCBI's set. Covers only versions a release
  annotated — not the older ones NCBI supplements the set with — and multiplies the store's size by the number of
  releases, when NCBI has already chosen the most recent alignment for each version.

- **Leave the historical set's current versions out of the shard**, since the release shard above it answers for them.
  Halves the shard, but makes it depend on being stacked under a particular release to be complete; a shard that holds
  the whole set is the set.

- **Keep a retired version's tags as history.** A tag is read as a recommendation; a reader would have to know the
  version's status to know whether it still held, which is what the status field is for.

## Open questions

- Whether the bundle should carry the successor Entrez names for a replaced version, so a caller could offer the current
  model's answer beside the retired one. The table has the column; the schema has no field yet.
- A store that stacks the historical set under GRCh38 but has nothing comparable for GRCh37 answers the same old name on
  one assembly and not the other. Whether the GRCh37 store should say so, rather than "not in the reference data", is
  open.
