#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
XS_PATH="${SCRIPT_DIR}/openmc_data/jeff311-hdf5/cross_sections.xml"
OUTPUT_DIR="${SCRIPT_DIR}/jeff3.1.1"

if [[ ! -r "${XS_PATH}" ]]; then
  echo "cross_sections.xml not found at ${XS_PATH}" >&2
  echo "Run ${SCRIPT_DIR}/convert_jeff311_ace.py first." >&2
  exit 1
fi

export OPENMC_CROSS_SECTIONS="${XS_PATH}"
export SPX_OUTPUT_DIR="${OUTPUT_DIR}"
export SPX_TEMPERATURE_METHOD="nearest"
export SPX_TEMPERATURE_TOLERANCE="1000.0"
export SPX_PLOT_PIXELS="2000"
mkdir -p "${SPX_OUTPUT_DIR}"
echo "OPENMC_CROSS_SECTIONS=${OPENMC_CROSS_SECTIONS}"
echo "SPX_OUTPUT_DIR=${SPX_OUTPUT_DIR}"
echo "SPX_TEMPERATURE_METHOD=${SPX_TEMPERATURE_METHOD}"
echo "SPX_TEMPERATURE_TOLERANCE=${SPX_TEMPERATURE_TOLERANCE}"
echo "SPX_PLOT_PIXELS=${SPX_PLOT_PIXELS}"
echo "You can now run: micromamba run -n openmc python ${SCRIPT_DIR}/run-spx.py"
