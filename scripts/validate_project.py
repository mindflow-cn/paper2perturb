#!/usr/bin/env python3
"""Validate Paper2Perturb skill structure and Python syntax."""

from __future__ import annotations

import compileall
import re
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SKILLS_ROOT = PROJECT_ROOT / "skills"
VALID_NAME = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
EXPECTED_KEYS = {"name", "description"}


def read_frontmatter(path: Path) -> dict[str, str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != "---":
        raise ValueError("missing opening YAML delimiter")

    try:
        end = lines.index("---", 1)
    except ValueError as exc:
        raise ValueError("missing closing YAML delimiter") from exc

    fields: dict[str, str] = {}
    for line in lines[1:end]:
        if not line.strip():
            continue
        if ":" not in line:
            raise ValueError(f"invalid frontmatter line: {line!r}")
        key, value = line.split(":", 1)
        fields[key.strip()] = value.strip().strip('"').strip("'")
    return fields


def validate_skill(skill_dir: Path) -> list[str]:
    errors: list[str] = []
    name = skill_dir.name
    skill_md = skill_dir / "SKILL.md"
    agent_yaml = skill_dir / "agents" / "openai.yaml"

    if not VALID_NAME.fullmatch(name):
        errors.append(f"{name}: directory name must be lowercase hyphen-case")
    if not skill_md.is_file():
        return errors + [f"{name}: missing SKILL.md"]

    try:
        fields = read_frontmatter(skill_md)
    except ValueError as exc:
        return errors + [f"{name}: {exc}"]

    if set(fields) != EXPECTED_KEYS:
        errors.append(
            f"{name}: frontmatter keys must be {sorted(EXPECTED_KEYS)}, "
            f"found {sorted(fields)}"
        )
    if fields.get("name") != name:
        errors.append(f"{name}: frontmatter name is {fields.get('name')!r}")
    if not fields.get("description"):
        errors.append(f"{name}: description is empty")

    if not agent_yaml.is_file():
        errors.append(f"{name}: missing agents/openai.yaml")
    elif f"${name}" not in agent_yaml.read_text(encoding="utf-8"):
        errors.append(f"{name}: default prompt must reference ${name}")

    for readme in skill_dir.rglob("README*"):
        errors.append(f"{name}: move per-skill documentation to the root README: {readme}")
    return errors


def main() -> int:
    skill_dirs = sorted(path for path in SKILLS_ROOT.iterdir() if path.is_dir())
    if not skill_dirs:
        print("ERROR: no skill directories found", file=sys.stderr)
        return 1

    errors = [error for skill in skill_dirs for error in validate_skill(skill)]
    if not compileall.compile_dir(PROJECT_ROOT / "skills", quiet=1):
        errors.append("one or more Python files failed to compile")
    if not compileall.compile_dir(PROJECT_ROOT / "scripts", quiet=1):
        errors.append("one or more project scripts failed to compile")

    if errors:
        print("Project validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    print(f"Validated {len(skill_dirs)} skills and all Python files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
