# weaver-data-provider development notes

## Product

See [`docs/PRODUCT.md`](docs/PRODUCT.md) for the product north star — what the package is, the load-bearing principles
(the publisher's placement, never a derived one; the files at rest are a contract; an answer or an error, never a
guess), and what is out of scope. Read it before proposing designs or plans. Shared terminology lives in
[`GLOSSARY.md`](GLOSSARY.md).

## Working norms

Operating directives for Claude (and any agent) in this repo; they counteract default model dispositions.

- **Resist the minimal-diff reflex.** Don't reach for the smallest change that hides the symptom (special-casing,
  papering over root causes). Aim for the correct fix at the right complexity level — not the smallest, not gold-plated.
- **Fail loudly and early.** Raise on a missing expected input or precondition; never fall back to a default/placeholder
  to limp along. A placeholder is an explicit caller input, never a code default.
- **Never instruct around a defect — fix the defect.** Don't write prose telling readers to work around broken code —
  "pass it as a string, the converter loses precision". Prose is untested, and callers who didn't read it stay broken.
- **Push back; don't just comply.** When a design, name, or approach seems worse — including a shortcut you're asked to
  take — say so with reasoning, unprompted. The author owns the final call.
- **Offer better alternatives with trade-offs.** When a materially better approach than the proposed one exists, present
  it and the trade-offs — don't just execute the ask.
- **Investigate before producing.** Read the code and verify constraints first. Don't treat a training-pattern
  convention as load-bearing unchecked; don't speculate about what you can read.
- **Explain non-obvious changes first.** For a change whose rationale isn't self-evident, give the why before showing or
  applying the diff.
- **Ask when unsure** rather than assume intent.
- **No intensifiers or emphasis filler.** Drop words and phrases that add emphasis but no information — "that's the
  key", "crucially", "importantly", "the key insight", "it's worth noting". State the point plainly. Applies to all
  prose: chat replies, PR/review comments, commit messages, and docs.

## Code style

@docs/style/general.md

Tests — what one asserts, what it may depend on, what its data may contain — follow
[`docs/style/writing-tests.md`](docs/style/writing-tests.md); it loads when a test file is touched.

## Schema

The protos under `proto/` are the source of truth for the files at rest; the committed Python stubs are generated. After
editing a proto: `buf lint`, `buf breaking --against .git#branch=main`, then
`uv run --group codegen python scripts/regen.py` and commit the stubs with the proto. Evolution is additive only — CI's
`schema-compat` job fails a renumber or type change. A change a reader cannot ignore bumps `FORMAT_VERSION`, so an old
reader refuses the new files rather than misreading them.

## Dependency direction

At runtime this package depends on hgvs-weaver and never the reverse ([`docs/PRODUCT.md`](docs/PRODUCT.md)). Don't
propose a change that makes weaver's runtime or CLI import this package; its test suite may.

## Docs

Two audiences, two registers:

- **Instruction files** are prompts and rules — `CLAUDE.md`, `.claude/rules/`, `.claude/skills/`: model-only, only what
  changes behavior, no maintainer notes, no harness mechanics (which rules load when, where files live). A token there
  is paid on every run that loads it; human-facing explanation belongs in `docs/` or code.
- **Everything under `docs/`** is written for a human first — a maintainer who has read
  [`docs/PRODUCT.md`](docs/PRODUCT.md) and [`GLOSSARY.md`](GLOSSARY.md) but not this area, and has to get the take-aways
  from one read on GitHub. Explain with the clarity and style of Martin Kleppmann — motivation before mechanism,
  specifics out of the argument's way. Detail that restates code — field lists, paths, env vars, test names — stays in
  the code and is linked, never transcribed. A model reads what a human reads. Design docs are the durable design
  record: one living doc per area under `docs/design/`, rewritten in place; no ADRs — rationale lives in the doc,
  chronology in git. The guide is [`docs/style/design-docs.md`](docs/style/design-docs.md); to write or rewrite one,
  load the `writing-design-docs` skill.

## Committing

- **Stage explicit paths**, not `git add -A` / `.`; explicit staging avoids sweeping in an untracked file.
- **Pre-commit runs lint/format/hygiene** (`.pre-commit-config.yaml`); CI runs the same hooks plus pytest. Ensure hooks
  are installed (`pre-commit install`) — if not, install or ask the author; never bypass with `--no-verify`.
- **Correct a pushed branch with a new commit on top**, not amend + force-push. PRs squash-merge, so `main` history
  stays linear regardless and intermediate fixups vanish on merge. Reserve force-push for rebasing a branch onto `main`.

## Worktrees

Worktrees go in `.claude/worktrees/` (gitignored), never `../` siblings.

- **New branch** → the Claude Code worktree command.
- **Existing branch** → `git worktree add .claude/worktrees/<name> <branch>` (the command only cuts fresh branches).

## CI and review

- **Adversarially review before opening a PR.** For any change with non-trivial code or logic, run adversarial review
  passes in subagents with fresh context — the reviewer sees only the diff, not the authoring conversation — and fix the
  findings autonomously; repeat until a pass surfaces only diminishing findings, then open the PR. Exempt: trivial
  changes, doc-only changes, resource/asset changes.
- **A PR description is written for the human reviewer**: what the change is and why, the take-aways, and where to look
  — the altitude of a design doc's Overview, shorter. The diff carries the detail; don't narrate it. Same style:
  [`docs/style/design-docs.md` § Style](docs/style/design-docs.md#style).
- **Pin third-party GitHub Actions to the latest stable release**: the moving major tag (`@v7`) where the action
  publishes one, else the exact latest version (`@v10.2.0`). Verify against the action's releases when adding or bumping
  one.
