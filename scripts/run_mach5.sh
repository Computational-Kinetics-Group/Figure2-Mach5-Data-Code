#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$root_dir/build" "$root_dir/runs" "$root_dir/output"

g++ -O3 -std=c++17 "$root_dir/src/order8_m5_validation.cpp" -o "$root_dir/build/order8_m5_validation"

# Four independent realizations; each produces 20 block estimates.  They are
# independent and can be run concurrently. Arguments after dt reproduce the
# publication protocol explicitly.
for seed in 1 2 3 4; do
  "$root_dir/build/order8_m5_validation" \
    5 "$seed" "$root_dir/runs/mach5_seed${seed}" 0.0003 \
    1.5 1.5 18000 300 12 20 \
    >"$root_dir/runs/mach5_seed${seed}.log" 2>&1 &
done
wait

python3 "$root_dir/scripts/aggregate_blocks.py" \
  "$root_dir"/runs/mach5_seed*_blocks.csv \
  -o "$root_dir/output/figure2_mach5_collision_production_regenerated.csv"

python3 "$root_dir/scripts/plot_figure2.py" \
  "$root_dir/output/figure2_mach5_collision_production_regenerated.csv" \
  -o "$root_dir/output/figure2_mach5_regenerated.png"
