# weaver-data-provider

Reference data for [hgvs-weaver](https://github.com/populationgenomics/hgvs-weaver): a builder that turns a publisher's
release into a local store, and the `DataProvider` that feeds weaver from it.

An HGVS engine is only as correct as what it is told about each transcript. A RefSeq transcript's sequence is not always
the genome's — about one record in forty differs — and where they differ by an insertion or deletion, where the
difference sits is a choice. NCBI publishes its choice: the alignment of every RefSeq transcript to the assembly. The
builder takes that alignment exon by exon rather than deriving one, so the positions weaver computes are the ones
ClinVar's names and VariantValidator's projections follow.

Every placement NCBI publishes is kept — on a chromosome, an alternate locus or a patch — so a gene the primary assembly
lacks is served where it is. An Ensembl transcript is placed at its exons, after a check that its record is the genome
spliced there. [`docs/design/placements.md`](docs/design/placements.md) has the reasons. A retired RefSeq version — the
one a variant name from years ago cites — is served from NCBI's historical set, stacked under the current release
([`docs/design/retired-versions.md`](docs/design/retired-versions.md)).

A store is one record per gene — every transcript, its protein, its alignments and the sequences they cite — in
[bagz](https://github.com/google-deepmind/bagz) files with sorted key and interval tables that are read into memory at
open. A lookup is a search over bytes in memory and one ranged read, from local disk or `gs://`.

```sh
weaver-data-build refseq ...     # one RefSeq release's bundles, as a shard
weaver-data-build historical ... # NCBI's set of retired RefSeq versions, as a shard to stack under a release's
weaver-data-build ensembl ...    # one Ensembl release's, placed at its exons and checked against the genome
weaver-data-build index ...      # the index and manifest over an ordered list of shards
weaver-data-build genome ...     # the assembly, cut into compressed blocks
```

The at-rest format is defined by the protos under `proto/`, and `buf breaking` gates every change: stores sit in other
projects' buckets, so the format only grows.

## Platforms

bagz publishes wheels for Linux x86_64 and macOS arm64; anywhere else, installing this package builds bagz from source.

## Development

`uv sync`, then `uv run pytest`. The Python stubs under `src/` are generated from `proto/` by
`uv run python scripts/regen.py` and committed.
