#!/usr/bin/env bash
# Fetch the large-tier instances from the GitHub Release.
#
# Usage: scripts/download_large.sh
#
# Files are downloaded into instances/ so all paths in grade.py and the
# README just work without any rewrite.

set -euo pipefail

REPO="${RELEASE_REPO:-ythuang0522/shortest-path-competition}"
TAG="${RELEASE_TAG:-v1.0}"

cd "$(dirname "$0")/.."
mkdir -p instances

files=(
  road_large.graph     road_large.queries
  grid_large.graph     grid_large.queries
  cluster_large.graph  cluster_large.queries
  social_large.graph   social_large.queries
)

echo "Downloading large-tier instances from ${REPO} @ ${TAG} ..."
for f in "${files[@]}"; do
  if [[ -f "instances/${f}" ]]; then
    echo "  [skip] instances/${f} already present"
    continue
  fi
  url="https://github.com/${REPO}/releases/download/${TAG}/${f}"
  echo "  [get ] ${f}"
  curl -fL --retry 3 -o "instances/${f}" "${url}"
done
echo "Done."
