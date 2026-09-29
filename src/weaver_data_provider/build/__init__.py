"""Build a store and a genome from a publisher's release, as local files.

Where the inputs come from and where the outputs go is the caller's: every function here reads local
paths and writes local paths, and records each input by role and by the digest of the bytes it read.
"""

from __future__ import annotations

import protovalidate
from google.protobuf import message as proto_message

import weaver_data_provider

BUILDER = f'weaver-data-provider {weaver_data_provider.__version__}'


class BuildError(RuntimeError):
    """A release this builder refuses to cut, rather than cut wrong."""


def validated(message: proto_message.Message, what: str) -> None:
    """Check a record against its proto's constraints before it is written.

    Raises:
        BuildError: Naming every violated field, the rule it broke and the value that broke it, since
            protovalidate's own message names only the message type.
    """
    try:
        protovalidate.validate(message)
    except protovalidate.ValidationError as error:
        details = []
        for violation in error.violations:
            field = violation.proto.field  # absent for a message-level rule
            path = '.'.join(e.field_name for e in field.elements) if field is not None else ''
            rule = f'{violation.proto.message} ({violation.proto.rule_id})'
            details.append(f'{path or "(message)"}: {rule}; got {violation.field_value!r}')
        raise BuildError(f'{what}: ' + '; '.join(details)) from error
