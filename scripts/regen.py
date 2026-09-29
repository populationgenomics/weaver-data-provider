"""Regenerate the committed Python stubs from proto/.

`buf export` materialises the protos with their buf.lock-pinned dependencies (buf/validate), then
grpcio-tools' protoc emits message classes and .pyi stubs into src/. buf/validate is emitted too: the
protovalidate wheels ship no Python stub for it, and the generated modules import it.

    uv run python scripts/regen.py
"""

from __future__ import annotations

import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / 'src'
GENERATED = (SRC / 'weaver_data_provider' / 'v1', SRC / 'buf')


def main() -> None:
    for path in GENERATED:
        if path.exists():
            shutil.rmtree(path)
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(['buf', 'export', str(ROOT), '-o', tmp], check=True)  # noqa: S603, S607
        protos = sorted(str(p.relative_to(tmp)) for p in pathlib.Path(tmp).rglob('*.proto'))
        subprocess.run(  # noqa: S603
            [sys.executable, '-m', 'grpc_tools.protoc', f'-I{tmp}', f'--python_out={SRC}', f'--pyi_out={SRC}', *protos],
            check=True,
        )
    for package in (SRC / 'weaver_data_provider' / 'v1', SRC / 'buf', SRC / 'buf' / 'validate'):
        (package / '__init__.py').touch()
    print('regenerated:', ', '.join(str(p.relative_to(ROOT)) for p in GENERATED))


if __name__ == '__main__':
    main()
