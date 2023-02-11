import numpy as np
from collections import Counter
import tensorflow as tf


def row_counter(my_array):
    list_of_tups = [tuple(el) for el in my_array]
    return Counter(list_of_tups)


def hist_labels(hist, CG_params):
    labels = []

    for enc in hist.keys():
        if CG_params['h_embed'] and CG_params['use_logits']: # and CG_params['use_STE'] is False:
            xs = np.array(enc)
        else:
            xs = 0.5*(np.array(enc) + 1)
            
        labels.append("".join(str(int(x)) for x in xs))
        #print("".join(str(int(x)) for x in xs))

    return labels


def genbitstrings(n):
    bitstrings = []

    for i in range(2**n):
        bs = np.binary_repr(i)
        bitstrings.append((n-len(bs))*'0' + bs)

    return bitstrings


def hamming2(s1, s2):
    """Calculate the Hamming distance between two bit strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def embedding_statistics(Hs, CG_params, hidden_dim=None):
    if hidden_dim is None:
        hidden_dim = CG_params['hidden_dim']

    if tf.is_tensor(Hs):
        Hs = Hs.numpy()

    hist = row_counter(Hs.astype(int))

    if CG_params['use_STE'] is False:
        cleaned_Hs = np.sign(Hs)
        hist = row_counter(cleaned_Hs)

    #print(hist)
    #print()

    labels = hist_labels(hist, CG_params)
    all_labels = genbitstrings(hidden_dim)

    counts = dict()

    for label in all_labels:
        if label in labels:

            handle = []
            for letter in label:
                if CG_params['h_embed'] and CG_params['use_logits']:
                    handle.append(int(letter))
                else:
                    handle.append(np.sign(int(letter)-0.5))

            counts[label] = hist[tuple(handle)]
        else:
            counts[label] = 0

    return counts, labels, all_labels, hist
