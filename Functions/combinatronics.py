# External
import itertools


def CartesianProduct(pots):
    cartesianProduct = []

    for j in itertools.product(*pots):
        cartesianProduct.append(j)

    return cartesianProduct
