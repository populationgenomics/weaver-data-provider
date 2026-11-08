"""The weaver data provider over a store and its genome: what an HGVS engine asks for, answered locally.

weaver's `VariantMapper` takes a `DataProvider`: transcript models with per-exon alignments, sequence
slices of transcripts, proteins and the genome, symbol and accession crosswalks, and refget
identities. One provider serves one assembly; a resolver working across builds holds one per build.
Every answer comes from the store's in-memory tables and one bundle read, or one indexed read of the
genome.
"""

from __future__ import annotations

import re
from collections.abc import Iterable

import weaver

from weaver_data_provider import genome as genome_mod
from weaver_data_provider import store as store_mod
from weaver_data_provider.v1 import bundle_pb2

# The reference publishers' namespaces, for an accession the store does not hold: it is still that kind, and
# weaver's error says it is not in the reference data rather than that the name is unknown.
_GENOMIC = ('NC_', 'NT_', 'NW_')
_TRANSCRIPT = ('NM_', 'NR_', 'XM_', 'XR_', 'ENST')
_PROTEIN = ('NP_', 'XP_', 'YP_', 'ENSP')
# UniProt's accession format, as its manual gives it, with an isoform's suffix or a sequence version: a bare
# accession in that shape is also the shape of gene symbols (H2AZP1) no store need hold.
_UNIPROT = re.compile(r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+|\.\d+)$')

# Whose records a gene's or a region's answer lists unless the caller says otherwise.
REFERENCE_PUBLISHERS = frozenset({bundle_pb2.PUBLISHER_REFSEQ, bundle_pb2.PUBLISHER_ENSEMBL})
_PUBLISHERS = frozenset(bundle_pb2.Publisher.values()) - {bundle_pb2.PUBLISHER_UNSPECIFIED}

# Whose accession a digest names first, when several publishers hold the residues; a protein no transcript
# encodes, UniProt's, after all of these. The order keeps the answer stable as shards are stacked.
_DIGEST_PREFERENCE = {bundle_pb2.PUBLISHER_REFSEQ: 0, bundle_pb2.PUBLISHER_ENSEMBL: 0, bundle_pb2.PUBLISHER_GENBANK: 1}
_UNRANKED = len(_DIGEST_PREFERENCE)


def _kind_name(kind: str | weaver.IdentifierType) -> str:
    return str(kind).lower()


def _sequence(bundle: bundle_pb2.GeneBundle, digest: str) -> bytes:
    for sequence in bundle.sequences:
        if sequence.digest == digest:
            return sequence.residues
    raise weaver.DataProviderError(f'{bundle.gene.symbol}: bundle cites sequence {digest} it does not carry')


class BundleProvider:
    """A `weaver.DataProvider` (and `Refget`, `TranscriptSearch`) over one assembly's store and genome.

    Every record the store holds is served when a name cites it. The transcripts and proteins listed for a
    gene, and the transcripts listed over a region, are by default RefSeq's and Ensembl's current ones: a
    GenBank record or a retired version would add a candidate the name did not cite, and weaver's readings
    say neither who published a record nor, without a catalogue, its status, so a caller could not take it
    out again. A caller wanting them listed says so when it makes the provider.
    """

    def __init__(
        self,
        store: store_mod.BundleStore,
        genome: genome_mod.Genome,
        *,
        listed_publishers: Iterable[bundle_pb2.Publisher] = REFERENCE_PUBLISHERS,
        include_retired: bool = False,
    ) -> None:
        """One assembly's provider.

        Args:
            store: The store, on the genome's assembly.
            genome: The genome.
            listed_publishers: Whose transcripts, and the proteins they encode, a gene's or a region's answer
                lists; a protein no transcript encodes is never listed. A record from any publisher is served
                when named.
            include_retired: List retired versions, superseded or suppressed, for a gene and over a region too.

        Raises:
            ValueError: If the store and the genome are on different assemblies, or `listed_publishers` is empty
                or names an undefined publisher.
        """
        if store.assembly != genome.assembly:
            name = bundle_pb2.Assembly.Name
            raise ValueError(f'store is on {name(store.assembly)}, genome on {name(genome.assembly)}')
        self._listed_publishers = frozenset(listed_publishers)
        if not self._listed_publishers or not self._listed_publishers <= _PUBLISHERS:
            raise ValueError(
                f'listed_publishers names no publisher, or one undefined: {sorted(self._listed_publishers)}'
            )
        self._store = store
        self._genome = genome
        self._include_retired = include_retired

    def _listed(self, transcript: bundle_pb2.Transcript) -> bool:
        """Whether a gene's or a region's answer lists this transcript."""
        current = transcript.status == bundle_pb2.TRANSCRIPT_STATUS_CURRENT
        return transcript.publisher in self._listed_publishers and (current or self._include_retired)

    def _gene_transcripts(self, symbol: str) -> list[tuple[bundle_pb2.GeneBundle, bundle_pb2.Transcript]]:
        return [(b, t) for b in self._store.by_symbol(symbol) for t in b.transcripts if self._listed(t)]

    # ---- what weaver reads ------------------------------------------------------------------------

    def _transcript(self, accession: str) -> tuple[bundle_pb2.GeneBundle, bundle_pb2.Transcript]:
        if '.' not in accession:
            raise weaver.DataProviderError(
                f'{accession}: a transcript model is per version; give the versioned accession'
            )
        base, _, version = accession.partition('.')
        # Records run across shards oldest release first, so the newest holding this version is last.
        for bundle in reversed(self._store.by_accession(accession)):
            for transcript in bundle.transcripts:
                if transcript.accession == base and str(transcript.version) == version:
                    return bundle, transcript
        raise weaver.DataProviderError(f'{accession}: not in the reference data')

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> weaver.TranscriptData:
        """The model on `reference_ac`, or on the first sequence it is placed on: a chromosome where there is one.

        A transcript placed on no sequence — a GenBank record, or a RefSeq one NCBI published no alignment
        for — is given without a placement whatever `reference_ac` is, so weaver reads it on its own sequence
        and refuses what needs a genome. A CDS end the publisher marks open is passed on as open, and weaver
        numbers no position from it.

        Raises:
            weaver.DataProviderError: If the transcript is not in the reference data, or is placed but not on
                `reference_ac`.
        """
        bundle, transcript = self._transcript(transcript_ac)
        coding = transcript.HasField('cds')
        model: weaver.TranscriptData = {
            'ac': transcript_ac,
            'gene': bundle.gene.symbol,
            'cds_start_index': transcript.cds.start_index if coding else None,
            'cds_end_index': transcript.cds.end_index_inclusive if coding else None,
            'cds_start_open': transcript.cds.start_open,
            'cds_end_open': transcript.cds.end_open,
        }
        if not transcript.alignments:
            return model
        alignments = [a for a in transcript.alignments if reference_ac is None or a.chromosome == reference_ac]
        if not alignments:
            raise weaver.DataProviderError(f'{transcript_ac}: no alignment on {reference_ac}')
        alignment = alignments[0]
        strand = 1 if alignment.strand == bundle_pb2.STRAND_PLUS else -1
        return {
            **model,
            'strand': strand,
            'reference_accession': alignment.chromosome,
            'exons': [
                {
                    'transcript_start': exon.transcript_start,
                    'transcript_end': exon.transcript_end,
                    'reference_start': exon.genome_start,
                    'reference_end': exon.genome_end_inclusive,
                    'alt_strand': strand,
                    'cigar': exon.cigar,
                }
                for exon in alignment.exons
            ],
        }

    def _is_protein(self, accession: str) -> bool:
        return self.get_identifier_type(accession) == weaver.IdentifierType.ProteinAccession

    def _residues(self, accession: str) -> bytes:
        if self._is_protein(accession):
            for bundle in reversed(self._store.by_protein(accession)):
                for protein in bundle.proteins:
                    if store_mod.protein_name(protein) == accession:
                        return _sequence(bundle, protein.sequence_digest)
            if '.' not in accession and self._store.holds_protein(accession):
                raise weaver.DataProviderError(
                    f'{accession}: a protein is held per version; give the versioned accession'
                )
            raise weaver.DataProviderError(f'{accession}: not in the reference data')
        bundle, transcript = self._transcript(accession)
        return _sequence(bundle, transcript.sequence_digest)

    def get_seq(self, ac: str, start: int, end: int | None, kind: str | weaver.IdentifierType) -> str:
        """Residues over the 0-based half-open range, clipped to the sequence at both ends."""
        if ac in self._genome:
            return self._genome.fetch(ac, start, end)
        if 'genomic' in _kind_name(kind):
            raise weaver.DataProviderError(f'{ac}: not a sequence of this genome')
        residues = self._residues(ac)
        start = max(start, 0)
        return residues[start:].decode('ascii') if end is None or end < 0 else residues[start:end].decode('ascii')

    def get_symbol_accessions(
        self, symbol: str, source_kind: str, target_kind: str
    ) -> list[tuple[weaver.IdentifierType, str]]:
        """`gene -> c` or `p -> c`: transcripts; `c -> p` or `gene -> p`: proteins; anything else is empty.

        A gene's proteins are those its listed transcripts encode. `p -> c` answers the transcripts recorded
        as encoding that version of the protein. A protein, like a transcript, is named by version: a bare
        accession, a version no release holds, or a protein no transcript encodes, a UniProt one, answers empty.
        """
        if source_kind == 'gene' and target_kind == 'c':
            transcripts = (f'{t.accession}.{t.version}' for _, t in self._gene_transcripts(symbol))
            return [(weaver.IdentifierType.TranscriptAccession, ac) for ac in dict.fromkeys(transcripts)]
        if source_kind == 'c' and target_kind == 'p':
            try:
                _, transcript = self._transcript(symbol)
            except weaver.DataProviderError:
                return []
            if not transcript.protein_accession:
                return []
            return [
                (weaver.IdentifierType.ProteinAccession, f'{transcript.protein_accession}.{transcript.protein_version}')
            ]
        if source_kind == 'p' and target_kind == 'c':
            if '.' not in symbol:
                return []
            base, _, version = symbol.partition('.')
            transcripts = (
                f'{t.accession}.{t.version}'
                for b in self._store.by_protein(symbol)
                for t in b.transcripts
                if t.protein_accession == base and str(t.protein_version) == version
            )
            return [(weaver.IdentifierType.TranscriptAccession, ac) for ac in dict.fromkeys(transcripts)]
        if source_kind == 'gene' and target_kind == 'p':
            proteins = (
                store_mod.protein_name(p)
                for b, t in self._gene_transcripts(symbol)
                for p in b.proteins
                if (p.transcript_accession, p.transcript_version) == (t.accession, t.version)
            )
            return [(weaver.IdentifierType.ProteinAccession, ac) for ac in dict.fromkeys(proteins)]
        return []

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        """What the store holds the identifier as; for one it does not hold, what its namespace says.

        UniProt's accession format, with an isoform's suffix or a version, is read last, after the gene symbols.
        """
        if identifier in self._genome or identifier.startswith(_GENOMIC):
            return weaver.IdentifierType.GenomicAccession
        if self._store.holds_transcript(identifier) or identifier.startswith(_TRANSCRIPT):
            return weaver.IdentifierType.TranscriptAccession
        if self._store.holds_protein(identifier) or identifier.startswith(_PROTEIN):
            return weaver.IdentifierType.ProteinAccession
        if self._store.by_symbol(identifier):
            return weaver.IdentifierType.GeneSymbol
        if _UNIPROT.match(identifier) is not None:
            return weaver.IdentifierType.ProteinAccession
        return weaver.IdentifierType.Unknown

    def get_transcripts_for_region(self, chrom: str, start: int, end: int) -> list[str]:
        """Versioned accessions of the listed transcripts placed over the 0-based closed interval, introns included.

        A retired version — one a historical shard holds as superseded or suppressed — is served when
        named, but not listed here unless the provider was made to: a region's answer is what the publisher
        recommends today.
        """
        found: dict[str, None] = {}
        for bundle in self._store.overlapping(chrom, start, end):
            for transcript in bundle.transcripts:
                if not self._listed(transcript):
                    continue
                if any(
                    a.chromosome == chrom
                    and min(e.genome_start for e in a.exons) <= end
                    and max(e.genome_end_inclusive for e in a.exons) >= start
                    for a in transcript.alignments
                ):
                    found[f'{transcript.accession}.{transcript.version}'] = None
        return list(found)

    # ---- refget ------------------------------------------------------------------------------------

    def get_refget_accession(self, ac: str) -> str | None:
        """The digest of the sequence the name serves: the newest release's, where several hold the name."""
        if ac in self._genome:
            return self._genome.digest(ac)
        try:
            if self._is_protein(ac):
                for bundle in reversed(self._store.by_protein(ac)):
                    for protein in bundle.proteins:
                        if store_mod.protein_name(protein) == ac:
                            return protein.sequence_digest
                return None
            _, transcript = self._transcript(ac)
        except weaver.DataProviderError:
            return None
        return transcript.sequence_digest

    def get_accession_for_refget(self, refget: str) -> str | None:
        """The accession carrying a digest: a chromosome, then a reference publisher's, then GenBank's, then UniProt's.

        A name with no version that a later release gives other residues, a reused UniProt isoform ID, no longer
        carries the digest and is not named.

        A `ga4gh:` prefix is accepted; a string that is not an `SQ.` digest names nothing, so answers None.
        """
        refget = refget.removeprefix('ga4gh:')
        if not refget.isascii() or not refget.startswith('SQ.'):
            return None
        chromosome = self._genome.name_for_digest(refget)
        if chromosome is not None:
            return chromosome
        try:
            holding = self._store.by_digest(refget)
        except ValueError:  # not a well-formed digest
            return None
        found: list[tuple[int, str]] = []
        for bundle in holding:
            publishers = {(t.accession, t.version): t.publisher for t in bundle.transcripts}
            for transcript in bundle.transcripts:
                if transcript.sequence_digest == refget:
                    rank = _DIGEST_PREFERENCE.get(transcript.publisher, _UNRANKED)
                    found.append((rank, f'{transcript.accession}.{transcript.version}'))
            for protein in bundle.proteins:
                if protein.sequence_digest == refget:
                    publisher = publishers.get((protein.transcript_accession, protein.transcript_version))
                    rank = _DIGEST_PREFERENCE.get(publisher, _UNRANKED) if publisher is not None else _UNRANKED
                    found.append((rank, store_mod.protein_name(protein)))
        # A name with no version, a UniProt isoform ID, may be one a later release gives other residues;
        # a versioned name names one sequence.
        for _, name in sorted(dict.fromkeys(found), key=lambda item: item[0]):
            if '.' in name or self.get_refget_accession(name) == refget:
                return name
        return None
