from buf.validate import validate_pb2 as _validate_pb2
from weaver_data_provider.v1 import bundle_pb2 as _bundle_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class KeyTable(_message.Message):
    __slots__ = ("kind", "keys", "key_offsets", "records", "record_offsets")
    KIND_FIELD_NUMBER: _ClassVar[int]
    KEYS_FIELD_NUMBER: _ClassVar[int]
    KEY_OFFSETS_FIELD_NUMBER: _ClassVar[int]
    RECORDS_FIELD_NUMBER: _ClassVar[int]
    RECORD_OFFSETS_FIELD_NUMBER: _ClassVar[int]
    kind: str
    keys: bytes
    key_offsets: bytes
    records: bytes
    record_offsets: bytes
    def __init__(self, kind: _Optional[str] = ..., keys: _Optional[bytes] = ..., key_offsets: _Optional[bytes] = ..., records: _Optional[bytes] = ..., record_offsets: _Optional[bytes] = ...) -> None: ...

class ChromosomeIntervals(_message.Message):
    __slots__ = ("assembly", "chromosome", "start", "end_inclusive", "record")
    ASSEMBLY_FIELD_NUMBER: _ClassVar[int]
    CHROMOSOME_FIELD_NUMBER: _ClassVar[int]
    START_FIELD_NUMBER: _ClassVar[int]
    END_INCLUSIVE_FIELD_NUMBER: _ClassVar[int]
    RECORD_FIELD_NUMBER: _ClassVar[int]
    assembly: _bundle_pb2.Assembly
    chromosome: str
    start: _containers.RepeatedScalarFieldContainer[int]
    end_inclusive: _containers.RepeatedScalarFieldContainer[int]
    record: _containers.RepeatedScalarFieldContainer[int]
    def __init__(self, assembly: _Optional[_Union[_bundle_pb2.Assembly, str]] = ..., chromosome: _Optional[str] = ..., start: _Optional[_Iterable[int]] = ..., end_inclusive: _Optional[_Iterable[int]] = ..., record: _Optional[_Iterable[int]] = ...) -> None: ...
