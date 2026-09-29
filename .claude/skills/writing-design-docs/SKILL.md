---
name: writing-design-docs
description: "Write or substantially rewrite a design doc under `docs/design/`. The procedure: what to read before drafting, what stays in the doc versus what moves to the code, the shape to draft in, the self-read, and the checks. Use when authoring a new design doc, rewriting an existing one, splitting or folding one, or bringing a doc up to the guide."
---

# Writing a design doc

The guidance is `docs/style/design-docs.md` — the reader, the style, where low-level detail goes, the default shape, the
policy. This is the procedure; it does not restate the guide.

## Read first

- The guide.
- What the doc's reader has already read: `docs/PRODUCT.md`, `GLOSSARY.md`.
- The doc as it stands, where one exists, and the Overview of each doc that will appear under `Related`.
- The code and contract files the doc describes — the proto, the module entry points, the tests. "The code states it" is
  verified there before a fact is left out, never assumed.

## Separate what stays from what restates code

Stays: the decisions, the named interfaces and what each promises, the consequences, the alternatives and why each was
rejected, open questions.

Goes: per-field paraphrase of a proto or schema, env-var names, paths beyond an entry point, function, class and test
names, error strings, constants.

For each code-owned fact the doc does not carry, note where it lives — the comment on the proto field, the docstring,
the test — or that it has no home yet because that code does not exist on this branch. Keep the list; it is part of the
report. No such fact goes into an inline comment beside the implementation (`docs/style/general.md`, Comments).

## Draft

The guide's default shape, deviating where the design reads better another way. A doc covering a surface, a flow, or
several interfaces uses the concrete forms the guide names: a mockup of the surface, a request diagram, a plain
statement of what is stored where, one subsection per interface in one consistent shape.

## Re-read as the reader

Read the draft as the maintainer the guide describes — has read `PRODUCT.md` and `GLOSSARY.md`, knows nothing about the
area, one read on GitHub. Can they state each decision and the reason for it afterwards? Where not, fix the passage, not
the reader.

## Checks

- `python3 tools/check_links.py` from the repo root — it scans all tracked Markdown, so read only the lines naming the
  doc.
- `pre-commit run mdformat --files <doc>`.
- After the push, check a mermaid block in the branch's file view on GitHub; the PR diff will not render it, so a
  reviewer opens the file.
- The header line carries `**Related:**`.

## Report

Where each code-owned fact lives, calling out the ones with no home yet so the author can push them onto the code's
documentation surface when it lands. Any question the writing raised for the author.
