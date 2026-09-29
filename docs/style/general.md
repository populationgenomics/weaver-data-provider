# Code style — language-agnostic

Principles independent of language; the language layer builds on it ([`python.md`](python.md)), and
[`writing-tests.md`](writing-tests.md) holds what a test asserts and may depend on. On conflict, the language doc wins
for its language. Code properties only — behavioral directives (push back, resist the minimal fix) live in
[`../../CLAUDE.md`](../../CLAUDE.md).

## Fail loud; never silently degrade

Raise on missing or malformed data; don't paper over it with `x or []`, `x or {}`, or a bare `if x:` that skips the real
case — that turns a missing input into a silent wrong answer. Validate and fail early.

"Not in the reference data" is a valid *answer* — a transcript version no release held
([`../PRODUCT.md`](../PRODUCT.md)) — not a license for code to swallow missing inputs: a lookup that finds nothing says
so; a store or genome that disagrees with its own manifest, or a publisher's file that breaks an expectation the builder
relies on, raises.

## Make the invalid state unrepresentable

Before writing a validation rule, ask why the value it rejects can be expressed at all. A field derived from another
cannot disagree with it; a type that admits only valid values needs no check; a record nested inside the thing it
belongs to cannot dangle. Each of those removes a rule, the silent failure the rule guarded against, and the test the
rule would have needed — for every caller at once, rather than for the paths a test happens to reach. A genome catalogue
that stored each sequence's first block number beside its length could disagree with itself; deriving the block number
from the lengths leaves nothing to disagree.

Fail loud ([above](#fail-loud-never-silently-degrade)) is the rule for the invalid input that can still arrive — from a
publisher's file, a bucket, a caller outside the type system. Between the two, prefer the shape that cannot hold the
error; validate what remains, and test that validation ([`writing-tests.md`](writing-tests.md)).

## Comments

Default to *no* comment. Add one only for a non-obvious *mechanism* or *constraint* a reader cannot recover from the
code — tersely, one line where possible. The *why* (why this shape was chosen) belongs in a design doc or the docstring,
not inline; never duplicate what a doc already states. A design-doc citation (`§N`, `see spike-infrastructure.md §3`)
*inside* a comment that also explains the design is the tell you are restating the doc — cut the explanation; a bare
one-line pointer is fine. No history narration ("removed X", "switched from Y"), no reference to a transient project
artifact ("this slice", "this spike", "this PR", "as X lands") — name the durable behavior, not the moment it arrived —
no commented-out code, no persuasion: write as if the current shape always existed. Self-check: a comment that stays
true after the code beneath it is rewritten is describing intent, not mechanism.

```
# Bad — rationale, persuasion, and the design doc already states this
conn = connect(dsn)  # one connection not a pool: pooling adds reconnect
# complexity we don't need yet, only pays off above N writers — the
# whole point of staying simple ...

# Good — one non-obvious fact; the why stays in the doc
conn = connect(dsn)  # single connection: the writer is single-threaded
```

```
# Bad — paraphrases the design doc and cites the section; intent, not mechanism
ipv4_enabled=True,  # Public IP, no authorized networks: reachable only through
# the connector (IAM-gated, TLS, ephemeral certs) — direct connections rejected.
# Private IP would need a VPC + serverless connector (spike-infrastructure.md §7).

# Good — the one non-obvious mechanism, terse
ipv4_enabled=True,  # empty authorizedNetworks ⇒ Cloud SQL refuses direct
# connections; the connector reaches it via an Admin-API ephemeral cert
```
