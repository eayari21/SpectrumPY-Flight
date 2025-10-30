#!/bin/bash

# || RUN SCIENCE: Run the science tool and backup the data

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}"
PYTHON_BIN="${PYTHON:-python}"

# Run science and optional backup routines through the packaged launcher
"$PYTHON_BIN" -m spectrumpy_flight.run_all
