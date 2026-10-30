"""Explicitly synthetic probability matrices for functional-defect tests."""
import numpy as np


def target_matrix(target):
    text = target.replace("+", "")
    matrix = np.eye(len(text))
    stack = []
    for i, symbol in enumerate(text):
        if symbol == "(":
            stack.append(i)
        elif symbol == ")":
            j = stack.pop()
            matrix[i, i] = matrix[j, j] = 0
            matrix[i, j] = matrix[j, i] = 1
    return matrix


def disrupt_pairs(matrix, indices, error):
    # Mix each perfect target pair with its unpaired state. Its two nucleotide
    # errors agree, preserving symmetry and row sums of probability one.
    pairs = {tuple(sorted((i, int(np.argmax(matrix[i]))))) for i in indices}
    for i, j in pairs:
        if i == j:
            raise AssertionError('Fixture expected paired target bases')
        matrix[i, j] = matrix[j, i] = 1-error
        matrix[i, i] = matrix[j, j] = error


def disrupt_unpaired(matrix, indices, error):
    # Introduce competing pair states inside this even-sized unpaired region.
    if len(indices) % 2:
        raise AssertionError('Fixture requires an even unpaired region')
    for i, j in zip(indices[::2], indices[1::2]):
        matrix[i, i] = matrix[j, j] = 1-error
        matrix[i, j] = matrix[j, i] = error
