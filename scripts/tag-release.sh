#!/bin/bash
set -e
set -o pipefail

ver=$1
if [ "$ver" == "" ]; then
	echo "Usage: $0 version" >&2
	exit 1
fi

nextnet_ver=$(cat extern/NEXTNet/.version)
case $ver in
	$nextnet_ver-*) ;;
	$nextnet_ver) ;;
	*)
		echo "Version $ver does not match $nextnet_ver" >&2
		exit 1
	;;
esac

if [ "$(git symbolic-ref --short HEAD)" != "master" ]; then
	echo "Currently checkout out branch must be 'master'" >&2
	exit 1
fi

# We ignore untracked files here
if [ $(git status -uno --porcelain | wc -l) != 0 ]; then
	echo "Working copy contains uncommitted changes" >&2
	exit 1
fi

if [ $(git ls-remote --tags origin v$ver | wc -l) != 0 ]; then
	echo "Version $ver already released (remote tag v$ver already exists)" >&2
	exit 1
fi

echo "Updating pyproject.toml" >&2
version = "0.4.0"
sed -i.bak 's/^version *=\(.*˜)$/version = "'"$ver"'"/' pyproject.toml
rm pyproject.toml.bak

echo "Comitting version bump to $ver" >&2
git add pyproject.toml
git commit -m "Incremented version to $ver"

echo "Tagging as v$ver and latest" >&2
git tag -f v$ver

echo "Pushing to origin" >&2
git push origin master v$ver

echo "Updating 'latest-release' on origin" >&2
git push -f origin v$ver:refs/tags/latest-release

