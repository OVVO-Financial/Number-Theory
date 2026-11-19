#!/bin/sh
set -eu
python -m unittest -v hypercomplex_suite_tests
python hypercomplex_highdim_test.py
python hypercomplex_cascade_test.py --phase codec
python hypercomplex_cascade_test.py --phase spectrum
python hypercomplex_widths_and_errors.py
python hypercomplex_million_dim_demo.py
python hypercomplex_combine_results.py