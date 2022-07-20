#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script contains the combinatronics functions
for creating the peptide libraries.
"""

import header as h


def cartesian_product(pots):
    """Calculates Cartesian Product to find all the
    unique combinations with which the peptide subunits can be
    combined."""

    cartesian_product_list = []

    # only creates combinations of amino acids listed in
    # their respective pots.
    for j in h.itertools.product(*pots):

        cartesian_product_list.append(j)

    return cartesian_product_list
