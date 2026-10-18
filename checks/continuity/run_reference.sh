#!/bin/sh
# Docker is the only runtime dependency on the host. Run from any directory.
set -eu
if [ "$#" -gt 1 ]; then
    echo "Usage: sh checks/continuity/run_reference.sh [NEW_OUTPUT_DIRECTORY]" >&2
    exit 2
fi
source_dir=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd -P)
output_dir=${1:-"$PWD/continuity-reference-output"}
# Require a new directory so earlier evidence cannot be overwritten.
mkdir -- "$output_dir"
output_dir=$(CDPATH= cd -- "$output_dir" && pwd -P)
image=continuity-hase-python3:pilot-01
docker build --tag "$image" "$source_dir/docker/continuity"
docker image inspect "$image" --format '{{.Id}}' > "$output_dir/image-id.txt"
docker run --rm --network none --read-only --cap-drop ALL \
    --security-opt no-new-privileges --cpus 2 --memory 1g --pids-limit 128 \
    --user "$(id -u):$(id -g)" --tmpfs /tmp:rw,nosuid,nodev,size=64m \
    --mount "type=bind,src=$source_dir,dst=/source,readonly" \
    --mount "type=bind,src=$output_dir,dst=/output" \
    --workdir /source "$image" sh -ec '
        python checks/continuity/hase_probe.py checks/continuity/fixture.json /output/observed.json
        python checks/continuity/compare_reference.py --observed /output/observed.json --output /output/comparison.json
    '
