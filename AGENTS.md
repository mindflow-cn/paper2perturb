# Repository Guidance

Paper2Perturb is the data curation agent for SimuCella. Keep this repository
independent from benchmark data and generated runtime outputs.

## Skill conventions

- Use lowercase hyphenated skill names and matching directory names.
- Keep only `name` and `description` in `SKILL.md` frontmatter.
- Put executable implementation in the skill's `scripts/` directory.
- Put detailed rules loaded on demand in `references/`.
- Do not add per-skill README files; update the root README files instead.
- Use repository-relative paths in documentation and scripts.
- Update `agents/openai.yaml` whenever a skill's purpose or name changes.

## Validation

Run these checks from the repository root after changing skills or scripts:

```bash
python3 scripts/validate_project.py
python3 skills/validate-metadata/scripts/validate_metadata.py --help
python3 skills/validate-h5ad/scripts/validate_h5ad.py --help
python3 skills/build-h5ad/scripts/convert.py --help
python3 skills/annotate-cell-types/scripts/test_annotate.py
```

Tests that download papers, models, or public datasets require network access
and should be run explicitly rather than as part of basic validation.
