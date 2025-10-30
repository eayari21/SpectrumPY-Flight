# Debugging the interactive clicking tools

The scatter-click launchers that live in `HDF_Explorer.py` and similar dashboards rely on matplotlib's picking events to dispatch the Quicklook viewer. When the callbacks stop responding, gathering a small amount of context up front lets us reproduce the issue locally and verify the fix quickly.

## Required information

Please collect the following details before opening an issue or sending a support request:

1. **Exact command used to launch the tool** – e.g., `python HDF_Explorer.py` from the repository root or a wrapper script. Include any flags or environment variables you set (such as `MPLBACKEND`).
2. **Dataset context** – the name of the HDF5 file that contains the selected event along with its full path relative to the repo (for example, `HDF5/example_run.h5`). Mention whether the file lives outside of the default `HDF5/` folder declared near the bottom of `HDF_Explorer.py`.【F:HDF_Explorer.py†L330-L474】
3. **Python environment** – interpreter version (`python --version`) and how the environment was created (conda env, venv, system Python). List any non-default matplotlib or Qt backends you installed.
4. **Operating system and display server** – macOS, Windows, or Linux, plus notable details like Wayland vs. X11 or remote/SSH sessions. These inputs matter because Qt and matplotlib vary their event-loop integrations by platform.
5. **GUI output** – copy the terminal log from the moment you launch the tool through the failed click. The Quicklook hand-off occurs inside the `HDF_Explorer.py` click callback; confirming that its log lines appear tells us whether the pick event fired.【F:HDF_Explorer.py†L478-L526】
6. **Reproduction steps** – describe the exact interaction sequence (e.g., "open the second scatter tab, hover over point 14, click once") and whether the point highlight, hover enlargement, or subprocess launch occurred. Mention whether you tried alternate files or datasets.

## Optional diagnostics

* Run `python -m pip list` (or `conda list`) to capture dependency versions.
* Toggle matplotlib's interactive backend by setting `MPLBACKEND=QtAgg` (or `MacOSX` on macOS) before launching the tool to see whether the issue is backend-specific.
* On Linux CI/remote machines, note whether you wrapped the launch command in `xvfb-run` or another virtual display helper.

Providing the above context lets us determine whether the click handler failed to register, the subprocess could not find `IDEX-quicklook.py`, or the environment blocked Qt from opening a new window. With a reproducible report, we can usually ship a patch in the next maintenance release.
