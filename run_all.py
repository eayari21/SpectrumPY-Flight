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
import os




def run_science_tool():
    try:
        science_tool_process = subprocess.Popen(["./science_tool.py"], shell=True)
        science_tool_process.wait()  # Wait for the subprocess to finish
    except KeyboardInterrupt:
        print("\nKeyboard interrupt received. Terminating subprocess.")
        science_tool_process.terminate()  # Terminate the subprocess if interrupted

def run_backup_script():
    subprocess.run(["./backup_science.py"], check=True, shell=True)

if __name__ == "__main__":
    run_science_tool()
    run_backup_script()