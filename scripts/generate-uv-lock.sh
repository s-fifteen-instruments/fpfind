#!/bin/sh
# Generate uv.lock file from poetry's pyproject.toml format
# Justin, 2025-12-01
#
# Alternative to use 'git stash; git restore pyproject.toml; git stash pop'.

# Go to project root
cd $(dirname $0)/..

# Preserve existing pyproject.toml
cp -p pyproject.toml pyproject.toml.bak

# Regenerate uv lockfile
rm -f uv.lock
uvx migrate-to-uv

# Restore pyproject.toml
mv pyproject.toml.bak pyproject.toml
