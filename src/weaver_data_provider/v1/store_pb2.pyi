from buf.validate import validate_pb2 as _validate_pb2
from weaver_data_provider.v1 import bundle_pb2 as _bundle_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class StoreManifest(_message.Message):
    __slots__ = ("format_version", "assembly", "shards", "builder", "index_files")
    FORMAT_VERSION_FIELD_NUMBER: _ClassVar[int]
    ASSEMBLY_FIELD_NUMBER: _ClassVar[int]
    SHARDS_FIELD_NUMBER: _ClassVar[int]
    BUILDER_FIELD_NUMBER: _ClassVar[int]
    INDEX_FILES_FIELD_NUMBER: _ClassVar[int]
    format_version: int
    assembly: _bundle_pb2.Assembly
    shards: _containers.RepeatedCompositeFieldContainer[Shard]
    builder: str
    index_files: _containers.RepeatedCompositeFieldContainer[IndexFile]
    def __init__(self, format_version: _Optional[int] = ..., assembly: _Optional[_Union[_bundle_pb2.Assembly, str]] = ..., shards: _Optional[_Iterable[_Union[Shard, _Mapping]]] = ..., builder: _Optional[str] = ..., index_files: _Optional[_Iterable[_Union[IndexFile, _Mapping]]] = ...) -> None: ...

class IndexFile(_message.Message):
    __slots__ = ("path", "digest", "digest_algorithm")
    PATH_FIELD_NUMBER: _ClassVar[int]
    DIGEST_FIELD_NUMBER: _ClassVar[int]
    DIGEST_ALGORITHM_FIELD_NUMBER: _ClassVar[int]
    path: str
    digest: str
    digest_algorithm: str
    def __init__(self, path: _Optional[str] = ..., digest: _Optional[str] = ..., digest_algorithm: _Optional[str] = ...) -> None: ...

class Shard(_message.Message):
    __slots__ = ("path", "release", "records", "digest", "digest_algorithm", "inputs")
    PATH_FIELD_NUMBER: _ClassVar[int]
    RELEASE_FIELD_NUMBER: _ClassVar[int]
    RECORDS_FIELD_NUMBER: _ClassVar[int]
    DIGEST_FIELD_NUMBER: _ClassVar[int]
    DIGEST_ALGORITHM_FIELD_NUMBER: _ClassVar[int]
    INPUTS_FIELD_NUMBER: _ClassVar[int]
    path: str
    release: str
    records: int
    digest: str
    digest_algorithm: str
    inputs: _containers.RepeatedCompositeFieldContainer[Input]
    def __init__(self, path: _Optional[str] = ..., release: _Optional[str] = ..., records: _Optional[int] = ..., digest: _Optional[str] = ..., digest_algorithm: _Optional[str] = ..., inputs: _Optional[_Iterable[_Union[Input, _Mapping]]] = ...) -> None: ...

class Input(_message.Message):
    __slots__ = ("role", "digest", "digest_algorithm")
    ROLE_FIELD_NUMBER: _ClassVar[int]
    DIGEST_FIELD_NUMBER: _ClassVar[int]
    DIGEST_ALGORITHM_FIELD_NUMBER: _ClassVar[int]
    role: str
    digest: str
    digest_algorithm: str
    def __init__(self, role: _Optional[str] = ..., digest: _Optional[str] = ..., digest_algorithm: _Optional[str] = ...) -> None: ...
