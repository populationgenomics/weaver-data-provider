# Retired transcript versions

**Related:** [`../PRODUCT.md`](../PRODUCT.md) — the principles this design applies; [`placements.md`](placements.md) —
how a transcript's placements and CDS are bundled, which a retired version shares;
[`../../GLOSSARY.md`](../../GLOSSARY.md) — annotation release, shard, store.

## Overview

A variant name cites a transcript version, and the name outlives the version's currency: ClinVar still carries
expressions on `NM_000059.3` years after NCBI moved to `.4`. An annotation release names only current versions, so a
store built from releases alone has no model for the names older variant records carry.

NCBI publishes a remedy: a historical set of every `NM_` and `NR_` version it has ever aligned to GRCh38 — about half of
them since replaced or suppressed — each with the alignment it last had, as the same kinds of file as a release. The
store takes that set as a shard of its own, stacked under the current release's shard, so a retired version is served
when named and the current release's model wins for a version both hold.

Each version carries NCBI's status for it — current, superseded or suppressed. The status is read not from the record,
which cannot state what happened after it was written, but from a table fetched from Entrez and given to the build. A
retired version carries no tags and no MANE partner, and a region's answer lists current versions only: a recommendation
is a property of today's release, not of a model that was once current.

## Background

A RefSeq transcript version is immutable: a change to the sequence is a new version, and the old one is *replaced*; a
record withdrawn is *suppressed*. The old version remains in GenBank and in Entrez, where its summary states its status
and, for a replaced one, its successor — which may be an unversioned accession, or another accession the record was
merged into. The record's own COMMENT states a replacement only when it happened before the record was last written, and
never states a suppression.

NCBI publishes the historical set (`RefSeq_historical_alignments`) under each GRCh38 annotation release's directory. It
is cumulative: every known version ever aligned to GRCh38, each with the alignment from the most recent release that
carried it, supplemented with older versions never annotated there, and filtered of alignments far from the gene's
expected location. It is anchored on one release, and the anchor can lag the release it sits beside: the set under
RS_2024_08 is anchored on RS_2023_03. It comes as an annotation GFF3 generated from the alignments, the alignment BAM
with the records' sequences, and the records as a GenBank flat file, but no sequence FASTA sets: the records carry each
transcript's sequence and, on its CDS feature, the protein id and translation. The GFF3 is generated from the alignments
rather than curated, and on about a quarter of its CDS rows it names another version's protein — `NP_056473.3` on
`NM_015658.3` and `NP_056473.2` on `NM_015658.4`, the two swapped — where a release's GFF3 never disagrees with the
record. A handful of its oldest records state their CDS in forms current records never use — a `complement` location on
an mRNA, a closed start read from its third base where `<` would say the start codon is missing. There is no such set
for GRCh37.

A store's manifest orders its shards, and the provider ([`provider.py`](../../src/weaver_data_provider/provider.py))
answers a named version from the last shard in that order that holds it. "Stacked under" means earlier in the order.

## Non-goals

- **No resolution of a retired version to its successor.** Entrez names the successor, and a caller wanting the current
  model's answer for an old name can ask for it; which version a study means is the study's decision.
- **No GRCh37 historical shard.** NCBI publishes no set for GRCh37; the only route is one shard per past release, which
  Alternatives weighs and rejects on size, and nothing yet justifies that cost there.
- **No status for an Ensembl transcript.** Ensembl archives releases rather than retiring versions; a retired Ensembl
  model is in an older release's shard, if a store stacks one.

## Design

### The historical set is a shard like any release's

The RefSeq builder reads the set's GFF3, BAM and GenBank files as it reads a release's, with differences the builder's
`historical` command ([`build/cli.py`](../../src/weaver_data_provider/build/cli.py)) expresses. The transcript and
protein sequences come from the GenBank records, since the set has no FASTA sets; a record with no sequence, or a coding
one whose CDS states no translation, fails the build. A status table is an input, and every version the annotation names
must have a row, so that each transcript's status is a fact the publisher stated rather than one the builder inferred.

The shard is stacked under the current release's shard. The store's rule that the last shard holding a version answers
for it then gives a version in both — one the release above still annotates — the current release's model, and a version
only the historical set holds, retired or unannotated, its last alignment. Nothing in the store changes for this: the
shard's immutability and the stacking order already carry it.

### The record names the protein

A transcript's protein is the one its GenBank record names, in the historical set and in a release alike; the
annotation's CDS rows name one too, and where the two differ the transcript is counted in the build's report. The record
is the publisher's statement about the transcript, and it is also what supplies the CDS and the translation, so the
protein accession, the protein sequence and the CDS the provider serves come from one source. The consequence is
visible: a `c.` variant on a retired version projects to a `p.` on the protein version that transcript version actually
encoded, which is the one ClinVar names beside it.

### A record the reader declines to interpret is left out, counted

The oldest records' legacy CDS forms are statements the reader does not interpret, and in a release such a record fails
the build, since a release is NCBI's current, consistent product and the form would be new. The historical set is an
archive, and five such records would refuse the whole set. In the historical set they are left out of the shard and
counted in the build's report with the reason, so the provider answers "not in the reference data" for them — a missing
answer, not a wrong one — and nothing about them is guessed.

### Status comes from Entrez, through a table

A version's status is a fact about the publisher's catalogue today. The record cannot state it, as Background says, and
the current annotation cannot either: a version absent from it may be live and unannotated (`NM_201563.5` on
RS_2024_08), and a version whose accession is absent may have been replaced rather than suppressed (`NM_000028.2`,
replaced by a `.3` that post-dates the release). Entrez states it, per version: live, replaced with a successor, or
suppressed; the three map onto `TranscriptStatus` in [`bundle.proto`](../../proto/weaver_data_provider/v1/bundle.proto).

The build reads local files only, so the Entrez answers are fetched beforehand by a command of their own,
`weaver-data-build status` ([`build/entrez.py`](../../src/weaver_data_provider/build/entrez.py)), into a table of one
row per version, and the table is an input to the build, recorded with the shard's other inputs by digest. A version
with no row, or a row with a status word that is not one of Entrez's, fails the build. The table is a snapshot: a
version live when the table was fetched and replaced since is current in the shard until the shard is cut again.

### A retired version carries no recommendation

The annotation's `tag` attribute and the MANE summary name what the publisher recommends: RefSeq Select, MANE Select. A
tag on a retired version in the set's GFF3 says what was recommended while that version was current, and the MANE
summary names current versions only. A retired version is bundled with neither, so a caller reading tags reads only
recommendations that hold. A version the table calls live but the release above does not annotate keeps the tags the
anchor release gave it; whether they still hold is not a fact the set states.

For the same reason the provider's region lookup lists current transcripts only. A retired version is served when named
— that is the shard's purpose — but a caller asking "which transcripts are here" is asking about today's annotation, and
the historical shard would otherwise answer with a decade of replaced models for every gene. The consequence for a
caller through weaver: it learns a version is retired only by not finding it in a region's answer, since weaver's
transcript model has no status field; the status is readable from the bundle.

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

- **Fail the build on a legacy CDS form**, as a release does. Five records would refuse 189,000; and **interpreting the
  form** — reading `complement` as the forward range, a closed start with `/codon_start=3` as open — would bundle a CDS
  the record did not state.

- **Keep a retired version's tags as history.** A tag is read as a recommendation; a reader would have to know the
  version's status to know whether it still held.

- **Take the protein from the annotation**, as the exon placements are. Right for a release, where the two sources
  agree, and wrong for 38,000 historical transcripts, each of which would then project its variants onto another
  version's protein.

- **List every version in a region's answer and let the caller filter by status.** The caller through weaver receives
  bare accessions and a transcript model with no status field, so it could not filter; it would meet a decade of
  replaced models for every gene with no way to tell them from the current ones.

## Open questions

- Whether the bundle should carry the successor Entrez names for a replaced version, so a caller could offer the current
  model's answer beside the retired one. The table carries the successor; the bundle does not.
- Whether a version's status should reach a caller through weaver, in the transcript model or beside it, so that a
  region's answer could list retired versions after all and let the caller decide.
- A store that stacks the historical set under GRCh38 but has nothing comparable for GRCh37 answers the same old name on
  one assembly and not the other. Whether the GRCh37 store should say so, rather than "not in the reference data", is
  open.
