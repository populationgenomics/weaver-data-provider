from buf.validate import validate_pb2 as _validate_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class Publisher(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    PUBLISHER_UNSPECIFIED: _ClassVar[Publisher]
    PUBLISHER_REFSEQ: _ClassVar[Publisher]
    PUBLISHER_ENSEMBL: _ClassVar[Publisher]

class TranscriptStatus(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    TRANSCRIPT_STATUS_UNSPECIFIED: _ClassVar[TranscriptStatus]
    TRANSCRIPT_STATUS_CURRENT: _ClassVar[TranscriptStatus]
    TRANSCRIPT_STATUS_SUPERSEDED: _ClassVar[TranscriptStatus]
    TRANSCRIPT_STATUS_SUPPRESSED: _ClassVar[TranscriptStatus]

class Tag(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    TAG_UNSPECIFIED: _ClassVar[Tag]
    TAG_MANE_SELECT: _ClassVar[Tag]
    TAG_MANE_PLUS_CLINICAL: _ClassVar[Tag]
    TAG_REFSEQ_SELECT: _ClassVar[Tag]
    TAG_ENSEMBL_CANONICAL: _ClassVar[Tag]

class Assembly(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    ASSEMBLY_UNSPECIFIED: _ClassVar[Assembly]
    ASSEMBLY_GRCH38: _ClassVar[Assembly]
    ASSEMBLY_GRCH37: _ClassVar[Assembly]

class Strand(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    STRAND_UNSPECIFIED: _ClassVar[Strand]
    STRAND_PLUS: _ClassVar[Strand]
    STRAND_MINUS: _ClassVar[Strand]

class AlignmentSource(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    ALIGNMENT_SOURCE_UNSPECIFIED: _ClassVar[AlignmentSource]
    ALIGNMENT_SOURCE_NCBI_BAM: _ClassVar[AlignmentSource]
    ALIGNMENT_SOURCE_MANE_PARTNER: _ClassVar[AlignmentSource]
    ALIGNMENT_SOURCE_COMPUTED: _ClassVar[AlignmentSource]

class Alphabet(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    ALPHABET_UNSPECIFIED: _ClassVar[Alphabet]
    ALPHABET_NUCLEOTIDE: _ClassVar[Alphabet]
    ALPHABET_PROTEIN: _ClassVar[Alphabet]
PUBLISHER_UNSPECIFIED: Publisher
PUBLISHER_REFSEQ: Publisher
PUBLISHER_ENSEMBL: Publisher
TRANSCRIPT_STATUS_UNSPECIFIED: TranscriptStatus
TRANSCRIPT_STATUS_CURRENT: TranscriptStatus
TRANSCRIPT_STATUS_SUPERSEDED: TranscriptStatus
TRANSCRIPT_STATUS_SUPPRESSED: TranscriptStatus
TAG_UNSPECIFIED: Tag
TAG_MANE_SELECT: Tag
TAG_MANE_PLUS_CLINICAL: Tag
TAG_REFSEQ_SELECT: Tag
TAG_ENSEMBL_CANONICAL: Tag
ASSEMBLY_UNSPECIFIED: Assembly
ASSEMBLY_GRCH38: Assembly
ASSEMBLY_GRCH37: Assembly
STRAND_UNSPECIFIED: Strand
STRAND_PLUS: Strand
STRAND_MINUS: Strand
ALIGNMENT_SOURCE_UNSPECIFIED: AlignmentSource
ALIGNMENT_SOURCE_NCBI_BAM: AlignmentSource
ALIGNMENT_SOURCE_MANE_PARTNER: AlignmentSource
ALIGNMENT_SOURCE_COMPUTED: AlignmentSource
ALPHABET_UNSPECIFIED: Alphabet
ALPHABET_NUCLEOTIDE: Alphabet
ALPHABET_PROTEIN: Alphabet

class GeneBundle(_message.Message):
    __slots__ = ("gene", "transcripts", "proteins", "sequences")
    GENE_FIELD_NUMBER: _ClassVar[int]
    TRANSCRIPTS_FIELD_NUMBER: _ClassVar[int]
    PROTEINS_FIELD_NUMBER: _ClassVar[int]
    SEQUENCES_FIELD_NUMBER: _ClassVar[int]
    gene: Gene
    transcripts: _containers.RepeatedCompositeFieldContainer[Transcript]
    proteins: _containers.RepeatedCompositeFieldContainer[Protein]
    sequences: _containers.RepeatedCompositeFieldContainer[Sequence]
    def __init__(self, gene: _Optional[_Union[Gene, _Mapping]] = ..., transcripts: _Optional[_Iterable[_Union[Transcript, _Mapping]]] = ..., proteins: _Optional[_Iterable[_Union[Protein, _Mapping]]] = ..., sequences: _Optional[_Iterable[_Union[Sequence, _Mapping]]] = ...) -> None: ...

class Gene(_message.Message):
    __slots__ = ("symbol", "hgnc_id", "ncbi_gene_id", "ensembl_gene_id", "previous_symbols", "alias_symbols", "name")
    SYMBOL_FIELD_NUMBER: _ClassVar[int]
    HGNC_ID_FIELD_NUMBER: _ClassVar[int]
    NCBI_GENE_ID_FIELD_NUMBER: _ClassVar[int]
    ENSEMBL_GENE_ID_FIELD_NUMBER: _ClassVar[int]
    PREVIOUS_SYMBOLS_FIELD_NUMBER: _ClassVar[int]
    ALIAS_SYMBOLS_FIELD_NUMBER: _ClassVar[int]
    NAME_FIELD_NUMBER: _ClassVar[int]
    symbol: str
    hgnc_id: str
    ncbi_gene_id: str
    ensembl_gene_id: str
    previous_symbols: _containers.RepeatedScalarFieldContainer[str]
    alias_symbols: _containers.RepeatedScalarFieldContainer[str]
    name: str
    def __init__(self, symbol: _Optional[str] = ..., hgnc_id: _Optional[str] = ..., ncbi_gene_id: _Optional[str] = ..., ensembl_gene_id: _Optional[str] = ..., previous_symbols: _Optional[_Iterable[str]] = ..., alias_symbols: _Optional[_Iterable[str]] = ..., name: _Optional[str] = ...) -> None: ...

class Transcript(_message.Message):
    __slots__ = ("accession", "version", "publisher", "status", "biotype", "tags", "sequence_digest", "cds", "protein_accession", "protein_version", "mane_partner", "alignments")
    ACCESSION_FIELD_NUMBER: _ClassVar[int]
    VERSION_FIELD_NUMBER: _ClassVar[int]
    PUBLISHER_FIELD_NUMBER: _ClassVar[int]
    STATUS_FIELD_NUMBER: _ClassVar[int]
    BIOTYPE_FIELD_NUMBER: _ClassVar[int]
    TAGS_FIELD_NUMBER: _ClassVar[int]
    SEQUENCE_DIGEST_FIELD_NUMBER: _ClassVar[int]
    CDS_FIELD_NUMBER: _ClassVar[int]
    PROTEIN_ACCESSION_FIELD_NUMBER: _ClassVar[int]
    PROTEIN_VERSION_FIELD_NUMBER: _ClassVar[int]
    MANE_PARTNER_FIELD_NUMBER: _ClassVar[int]
    ALIGNMENTS_FIELD_NUMBER: _ClassVar[int]
    accession: str
    version: int
    publisher: Publisher
    status: TranscriptStatus
    biotype: str
    tags: _containers.RepeatedScalarFieldContainer[Tag]
    sequence_digest: str
    cds: Cds
    protein_accession: str
    protein_version: int
    mane_partner: str
    alignments: _containers.RepeatedCompositeFieldContainer[Alignment]
    def __init__(self, accession: _Optional[str] = ..., version: _Optional[int] = ..., publisher: _Optional[_Union[Publisher, str]] = ..., status: _Optional[_Union[TranscriptStatus, str]] = ..., biotype: _Optional[str] = ..., tags: _Optional[_Iterable[_Union[Tag, str]]] = ..., sequence_digest: _Optional[str] = ..., cds: _Optional[_Union[Cds, _Mapping]] = ..., protein_accession: _Optional[str] = ..., protein_version: _Optional[int] = ..., mane_partner: _Optional[str] = ..., alignments: _Optional[_Iterable[_Union[Alignment, _Mapping]]] = ...) -> None: ...

class Cds(_message.Message):
    __slots__ = ("start_index", "end_index_inclusive")
    START_INDEX_FIELD_NUMBER: _ClassVar[int]
    END_INDEX_INCLUSIVE_FIELD_NUMBER: _ClassVar[int]
    start_index: int
    end_index_inclusive: int
    def __init__(self, start_index: _Optional[int] = ..., end_index_inclusive: _Optional[int] = ...) -> None: ...

class Alignment(_message.Message):
    __slots__ = ("assembly", "chromosome", "strand", "source", "exons")
    ASSEMBLY_FIELD_NUMBER: _ClassVar[int]
    CHROMOSOME_FIELD_NUMBER: _ClassVar[int]
    STRAND_FIELD_NUMBER: _ClassVar[int]
    SOURCE_FIELD_NUMBER: _ClassVar[int]
    EXONS_FIELD_NUMBER: _ClassVar[int]
    assembly: Assembly
    chromosome: str
    strand: Strand
    source: AlignmentSource
    exons: _containers.RepeatedCompositeFieldContainer[Exon]
    def __init__(self, assembly: _Optional[_Union[Assembly, str]] = ..., chromosome: _Optional[str] = ..., strand: _Optional[_Union[Strand, str]] = ..., source: _Optional[_Union[AlignmentSource, str]] = ..., exons: _Optional[_Iterable[_Union[Exon, _Mapping]]] = ...) -> None: ...

class Exon(_message.Message):
    __slots__ = ("transcript_start", "transcript_end", "genome_start", "genome_end_inclusive", "cigar")
    TRANSCRIPT_START_FIELD_NUMBER: _ClassVar[int]
    TRANSCRIPT_END_FIELD_NUMBER: _ClassVar[int]
    GENOME_START_FIELD_NUMBER: _ClassVar[int]
    GENOME_END_INCLUSIVE_FIELD_NUMBER: _ClassVar[int]
    CIGAR_FIELD_NUMBER: _ClassVar[int]
    transcript_start: int
    transcript_end: int
    genome_start: int
    genome_end_inclusive: int
    cigar: str
    def __init__(self, transcript_start: _Optional[int] = ..., transcript_end: _Optional[int] = ..., genome_start: _Optional[int] = ..., genome_end_inclusive: _Optional[int] = ..., cigar: _Optional[str] = ...) -> None: ...

class Protein(_message.Message):
    __slots__ = ("accession", "version", "sequence_digest", "transcript_accession", "transcript_version")
    ACCESSION_FIELD_NUMBER: _ClassVar[int]
    VERSION_FIELD_NUMBER: _ClassVar[int]
    SEQUENCE_DIGEST_FIELD_NUMBER: _ClassVar[int]
    TRANSCRIPT_ACCESSION_FIELD_NUMBER: _ClassVar[int]
    TRANSCRIPT_VERSION_FIELD_NUMBER: _ClassVar[int]
    accession: str
    version: int
    sequence_digest: str
    transcript_accession: str
    transcript_version: int
    def __init__(self, accession: _Optional[str] = ..., version: _Optional[int] = ..., sequence_digest: _Optional[str] = ..., transcript_accession: _Optional[str] = ..., transcript_version: _Optional[int] = ...) -> None: ...

class Sequence(_message.Message):
    __slots__ = ("digest", "alphabet", "residues")
    DIGEST_FIELD_NUMBER: _ClassVar[int]
    ALPHABET_FIELD_NUMBER: _ClassVar[int]
    RESIDUES_FIELD_NUMBER: _ClassVar[int]
    digest: str
    alphabet: Alphabet
    residues: bytes
    def __init__(self, digest: _Optional[str] = ..., alphabet: _Optional[_Union[Alphabet, str]] = ..., residues: _Optional[bytes] = ...) -> None: ...
