#!/usr/bin/env bash
# Publishes the Doxygen HTML of a version in the gh-pages branch (GitHub Pages: Deploy from a branch,
# gh-pages, / (root)), used by .github/workflows/docs.yml:
#
#   publish-docs.sh <html directory> <version: master or vX.Y.Z> <version-selector.js>
#
# The site keeps a folder per version (the one given is replaced) and, at its root, the version selector,
# versions.json (master, then the releases from the newest) and index.html, opening the newest release (master
# while there is none). The branch is a single commit, replaced at each publication (the history of the
# documentation is not kept); a concurrent publication (another version) makes the push fail: it is retried
# on the new branch. Needs GITHUB_TOKEN (contents: write), GITHUB_REPOSITORY and GITHUB_SHA; DOCS_REMOTE sets
# another repository (e.g. a local bare one to try the script).
set -euo pipefail

html=$(realpath "$1")
version=$2
selector=$(realpath "$3")
[[ "$version" =~ ^(master|v[0-9]+\.[0-9]+\.[0-9]+)$ ]] || {
    echo "Invalid version '$version': master or vX.Y.Z"
    exit 1
}

remote="${DOCS_REMOTE:-https://x-access-token:${GITHUB_TOKEN}@github.com/${GITHUB_REPOSITORY}.git}"
site=$(mktemp -d)

for attempt in 1 2 3 4 5; do
    rm -rf "$site" && mkdir -p "$site"
    # The current site (none at the first publication)
    old=$(git ls-remote "$remote" refs/heads/gh-pages | cut -f1)
    if [[ -n "$old" ]]; then
        git clone -q --depth 1 --branch gh-pages "$remote" "$site"
    else
        git -C "$site" init -q
    fi
    git -C "$site" checkout -q --orphan publication

    rm -rf "${site:?}/$version"
    cp -r "$html" "$site/$version"
    cp "$selector" "$site/version-selector.js"
    # Doxygen file names start with _: no Jekyll processing
    touch "$site/.nojekyll"

    versions=$(
        cd "$site"
        [[ -d master ]] && echo master
        find . -maxdepth 1 -type d -name 'v*' -printf '%f\n' | sort -rV
    )
    newest=$(grep -m1 '^v' <<<"$versions" || echo master)
    sed 's/.*/"&"/' <<<"$versions" | paste -sd, | sed 's/.*/[&]/' >"$site/versions.json"
    cat >"$site/index.html" <<EOF
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>REST documentation</title>
<meta http-equiv="refresh" content="0; url=${newest}/index.html">
</head>
<body><p><a href="${newest}/index.html">REST documentation ${newest}</a></p></body>
</html>
EOF

    git -C "$site" add -A
    git -C "$site" -c user.name="github-actions[bot]" \
        -c user.email="41898282+github-actions[bot]@users.noreply.github.com" \
        commit -q -m "Documentation ${version} (${GITHUB_SHA:0:7})"
    if git -C "$site" push -q --force-with-lease="gh-pages:${old}" "$remote" HEAD:refs/heads/gh-pages; then
        echo "Documentation ${version} published; versions: $(paste -sd' ' <<<"$versions")"
        exit 0
    fi
    echo "The gh-pages branch changed (attempt ${attempt}), retrying"
    sleep $((attempt * 10))
done
echo "The documentation could not be published"
exit 1
