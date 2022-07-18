#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script contains the combinatronics functions for creating the peptide libraries.
"""

from scripts import header as h


def CartesianProduct(pots):
    """Calculates Cartesian Product to find all the unique combinations with which the peptide subunits can be
    combined."""

    cartesian_product = []

    # only creates combinations of amino acids listed in their respective pots.
    for j in h.itertools.product(*pots):

        cartesian_product.append(j)

    return cartesian_product