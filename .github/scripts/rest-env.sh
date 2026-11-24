#!/usr/bin/env bash
# The environment of this commit (CONTRIBUTING.md, Continuous integration): the tag of the rest-docker
# repository in .github/rest-env, for its image ghcr.io/<owner>/rest-env:<tag> and its rest/Dockerfile.
# Written as the outputs ref and image of the step; the repository variable RESTV2_IMAGE (environment)
# replaces the image.
set -euo pipefail

ref=$(tr -d '[:space:]' < .github/rest-env)
if [[ ! "$ref" =~ ^v[0-9]+(\.[0-9]+)*$ ]]; then
  echo "The tag of rest-docker in .github/rest-env is not vN (or vN.M...): '$ref'"
  exit 1
fi
image=${RESTV2_IMAGE:-ghcr.io/${GITHUB_REPOSITORY_OWNER,,}/rest-env:${ref}}
echo "rest-docker $ref, image $image"
{
  echo "ref=$ref"
  echo "image=$image"
} >> "$GITHUB_OUTPUT"
