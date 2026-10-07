#!/usr/bin/env bash
# Run the three-body problem once for each integrator, then open all four
# results in the visualizer at the same time.
#
# Usage:
#   ./scripts/compare_threebodies.sh [t-end] [fps]
#
# t-end is the simulation duration passed to threebodies (default 2000).
# fps is the visualizer playback rate (default 60).

set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
threebodies="$root/build/apps/simulators/threebodies"
visualizer="$root/build/apps/visualizer/visualizer"
outdir="$root/analysis/results"

t_end="${1:-100000}"
fps="${2:-60}"

integrators=(
    velocity-verlet
)

if [[ ! -x "$threebodies" || ! -x "$visualizer" ]]; then
    echo "Build the simulators and visualizer first:" >&2
    echo "  cmake -B build && cmake --build build" >&2
    exit 1
fi

mkdir -p "$outdir"

json_files=()
for integrator in "${integrators[@]}"; do
    output="$outdir/threebodies_${integrator}.json"
    echo "simulating $integrator -> $output"
    "$threebodies" "$t_end" --integrator "$integrator" --output "$output"
    json_files+=("$output")
done

'''
pids=()
for output in "${json_files[@]}"; do
    echo "visualizing $output"
    "$visualizer" "$output" --fps "$fps" --trajectory-path-length 20000 &
    pids+=("$!")
done

wait "${pids[@]}"'''
