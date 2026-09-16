Policy on AI Agents
===================

This policy is experimental and is subject to change at the whims of the
maintainers.

Use of AI agents is permitted, but must be declared explicitly. Any submitted
contribution made with the assistance of an AI agent must have been carefully
checked by a qualified human. All contributions will additionally need to be
reviewed by one of the maintainers before merge, which we are happy to do, but
we do not have unlimited time to spend on low-effort, automated submissions.

Declare AI use in the commit message (a `Co-Authored-By:` trailer naming the
tool and version is appropriate) and in the pull request description. Reviewers
should be able to see at a glance that AI was involved and which tool was used.

"Carefully checked by a qualified human" means the human can explain why each
change is correct---not merely that it compiles and passes tests.
Rubber-stamping AI output without understanding it is not acceptable.

All contributions should pass both the `release` and `sanitize` builds and their
tests before submission, and `./run-tests.bash`, which is not registered with
`ctest` and so is run by neither the presets nor CI. Install
[VeriPB](https://gitlab.com/MIAOresearch/software/VeriPB) first: the
proof-verification tests are only registered when `veripb` is on the `PATH`, so
a green `ctest` without it has checked considerably less than it appears to.
See `README.md` for the build commands and `dev_docs/` for the design notes.

Licensing
=========

The solver is released under the terms in `LICENCE`, which is the MIT licence;
that file is the authoritative statement. Unless you explicitly state otherwise,
any contribution you intentionally submit for inclusion in this work shall be
licensed as above, without any additional terms or conditions. There is no
copyright assignment: you keep the copyright in what you write.

If a contribution contains code you did not write yourself, say where it came
from and under what terms in the pull request, so that we can check it can be
released before merging.

Code Formatting
===============

All C++ source is formatted with
[clang-format](https://clang.llvm.org/docs/ClangFormat.html) using the
`.clang-format` in the repository root. Format the whole tree before submitting:

```shell
git ls-files '*.cc' '*.hh' | xargs clang-format -i
```

Use **clang-format 21**: output can differ between major releases (18 and 21
disagree about the continuation indent of a multi-line `for` header, for one),
and CI pins 21.1.8, so a different major version will report changes that have
nothing to do with your contribution. The PyPI wheel is the easiest way to get
exactly it:

```shell
pip install clang-format==21.1.8
```

The `clang-format` workflow runs `clang-format --dry-run --Werror` over the same
files on every push and pull request, so an unformatted contribution fails the
check.

`clang-format` settles only the mechanical questions, and two conventions it
cannot settle are worth knowing:

- The empty `//` markers in the cxxopts option blocks in `src/` are
  load-bearing. Each block is a single chained expression, and with
  `ColumnLimit: 0` clang-format would pack it onto one line hundreds of
  characters long; a trailing comment cannot be joined with what follows it, so
  the marker pins the break. Keep them when editing those blocks, and add one to
  any new option line that a chained call follows. `.clang-format` carries a
  comment saying so.
- `AlignTrailingComments` is off deliberately. With no column limit, aligning
  comments across a block of long lines pushes them a very long way right.

Branches and History
====================

Work on a branch and open a pull request; do not commit to `main`.

A pull request's MERGED status is **not** a reliable signal that its content is
in `main`, and stacked pull requests are how that goes wrong. A PR based on
another branch rather than on `main` is merged *into that branch*, and GitHub
marks it MERGED at that point --- whether or not the parent ever carries it the
rest of the way. To decide whether a branch's content has landed, ask the
content rather than the PR status:

```shell
git merge-base --is-ancestor origin/<branch> origin/main   # exit 0: it is in
git cherry origin/main origin/<branch>                     # no '+' lines: it is in
```

Before deleting a remote branch, check that no open pull request is based on it
(`gh pr list --base <branch> --state open`). Deleting a branch that is an open
PR's base **closes that PR**; GitHub does not retarget it, and the close is not
reversible by recreating the branch.

Merge pull requests with a merge commit, which is what the history here uses.

Developer Documentation
=======================

Architectural notes on individual subsystems live in `dev_docs/`:
[`architecture.md`](dev_docs/architecture.md) for the layer map,
[`file-formats.md`](dev_docs/file-formats.md) for graph input and `InputGraph`,
[`proof-logging.md`](dev_docs/proof-logging.md) for the VeriPB machinery and
which option combinations support it, and
[`preprocessor-refactor.md`](dev_docs/preprocessor-refactor.md) for where the
pre-search code is going.

AI agents in particular should read the relevant document before making
non-trivial changes to a subsystem --- this is the most efficient way to absorb
the design decisions and conventions that are not obvious from the code alone.
