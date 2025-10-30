#!/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-

title = """

//  ================================================================================
//  ||                                                                            ||
//  ||              run_all                                                       ||
//  ||              ------------------------------------------------------        ||
//  ||                                  R U N  A L L                              ||
//  ||              ------------------------------------------------------        ||
//  ||                                                                            ||
//  ||                __author__      = Ethan Ayari                               ||
//  ||                IMPACT/LASP, CU Boulder                                     ||
//  ||                                                                            ||
//  ||                For: IDEX Flight Model Integration and Test, L0-L1A         ||
//  ||                                                                            ||
//  ||                2023                                                        ||
//  ||                                                                            ||
//  ||                                                                            ||
//  ||                Works with Python 3.10.4                                    ||
//  ||                                                                            ||
//  ================================================================================


run_all: A script that calls the science tool and the backup script.

"""

print(title)

# %%
# || Python libraries
import subprocess
import sys
from importlib import util


def run_science_tool():
    try:
        subprocess.run(
            [sys.executable, "-m", "spectrumpy_flight.science_tool"],
            check=True,
        )
    except KeyboardInterrupt:
        print("\nKeyboard interrupt received. Terminating subprocess.")


def run_backup_script():
    module_name = "spectrumpy_flight.backup_science"
    if util.find_spec(module_name) is None:
        print("Backup script not found; skipping backup stage.")
        return
    subprocess.run([sys.executable, "-m", module_name], check=True)


if __name__ == "__main__":
    run_science_tool()
    run_backup_script()
