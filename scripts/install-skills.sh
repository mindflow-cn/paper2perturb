#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET_ROOT="${1:-$(pwd)}"

if [[ ! -d "$TARGET_ROOT" ]]; then
    echo "ERROR: target directory does not exist: $TARGET_ROOT" >&2
    exit 1
fi

for client_dir in .agents/skills .claude/skills; do
    install_dir="$TARGET_ROOT/$client_dir"
    mkdir -p "$install_dir"
    for skill_dir in "$PROJECT_ROOT"/skills/*; do
        skill_name="$(basename "$skill_dir")"
        ln -sfn "$skill_dir" "$install_dir/$skill_name"
        echo "Installed $skill_name -> $install_dir/$skill_name"
    done
done
