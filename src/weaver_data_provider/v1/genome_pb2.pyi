from buf.validate import validate_pb2 as _validate_pb2
from weaver_data_provider.v1 import bundle_pb2 as _bundle_pb2
from weaver_data_provider.v1 import store_pb2 as _store_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class SequenceCatalogue(_message.Message):
    __slots__ = ("format_version", "assembly", "sequences", "source", "builder", "block_size", "digest", "digest_algorithm")
    FORMAT_VERSION_FIELD_NUMBER: _ClassVar[int]
    ASSEMBLY_FIELD_NUMBER: _ClassVar[int]
    SEQUENCES_FIELD_NUMBER: _ClassVar[int]
    SOURCE_FIELD_NUMBER: _ClassVar[int]
    BUILDER_FIELD_NUMBER: _ClassVar[int]
    BLOCK_SIZE_FIELD_NUMBER: _ClassVar[int]
    DIGEST_FIELD_NUMBER: _ClassVar[int]
    DIGEST_ALGORITHM_FIELD_NUMBER: _ClassVar[int]
    format_version: int
    assembly: _bundle_pb2.Assembly
    sequences: _containers.RepeatedCompositeFieldContainer[SequenceEntry]
    source: _store_pb2.Input
    builder: str
    block_size: int
    digest: str
    digest_algorithm: str
    def __init__(self, format_version: _Optional[int] = ..., assembly: _Optional[_Union[_bundle_pb2.Assembly, str]] = ..., sequences: _Optional[_Iterable[_Union[SequenceEntry, _Mapping]]] = ..., source: _Optional[_Union[_store_pb2.Input, _Mapping]] = ..., builder: _Optional[str] = ..., block_size: _Optional[int] = ..., digest: _Optional[str] = ..., digest_algorithm: _Optional[str] = ...) -> None: ...

class SequenceEntry(_message.Message):
    __slots__ = ("name", "length", "digest", "description")
    NAME_FIELD_NUMBER: _ClassVar[int]
    LENGTH_FIELD_NUMBER: _ClassVar[int]
    DIGEST_FIELD_NUMBER: _ClassVar[int]
    DESCRIPTION_FIELD_NUMBER: _ClassVar[int]
    name: str
    length: int
    digest: str
    description: str
    def __init__(self, name: _Optional[str] = ..., length: _Optional[int] = ..., digest: _Optional[str] = ..., description: _Optional[str] = ...) -> None: ...
