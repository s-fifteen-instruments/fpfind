#!/bin/sh
# Install pre-push hook to disable 'wip' commits
# Justin, 2025-12-02

# Go to project root
cd $(dirname $0)/..

# Copy pre-push hook into .git
cp -p scripts/hooks/pre-push .git/hooks/pre-push
