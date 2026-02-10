# SUDA FSW Science Packets White Paper (LASP Doc 168606)

This repository uses **LASP Document 168606, "SUDA FSW Science Packets White Paper" (Rev E)**
as the authority for packed TXHDR words and unpacked metadata derivations.

The parser updates in this branch annotate HDF5 metadata fields with this authority and
attach per-field descriptions for TXHDR packed and unpacked products.

If additional packet-table details are transcribed from the white paper, extend
`lookup/txhdr_descriptions.py`.
