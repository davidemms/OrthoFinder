#!usr/bin/python
import numpy as np

def wpgma(M, names = None, q_min=False):
    q_debug = False
    M = 0.5*( M + M.transpose())
    if q_min:
        np.fill_diagonal(M, 9e99)
        M = -M
    if names is None:
        forest = list(range(M.shape[0]))
    else:
        forest = names
    while len(forest) > 1:
        if q_debug: print("")
        i,j = get_max(M)
        if q_debug: print(i,j)
        M = combine(M, i, j)
        if q_debug: print(M)
        # remove j, combine (i,j) in i
        if i > j:
            swap = i
            i = j
            j = swap
        forest = forest[:i] + [(forest[i], forest[j]), ] + forest[i+1:j] + forest[j+1:]
        if q_debug: print(forest)
    return forest[0]

def get_max(M):
    return np.unravel_index(np.argmax(M, axis=None), M.shape)

def combine(M, i, j):
    r = 0.5*(M[i,:] + M[j,:])
    c = 0.5*(M[:,i] + M[:,j])
    r = r[0,:]
    M[i,:] = r
    M[:,i] = c
    M[i,i] = 0.
    M = np.delete(M, j, 0)
    M = np.delete(M, j, 1)
    return M

if __name__ == "__main__":
    # test (wikipedia)
    # Distance matrix
    M = np.matrix("0 17 21 31 23; 17 0 30 34 21; 21 30 0 28 39; 31 34 28 0 43; 23 21 39 43 0")
    # print(M)
    # x = wpgma(M, q_min=True)
    # assert((((0, 1), 4), (2, 3)) == x)

    # Similarity matrix
    M = np.matrix([[0.        , 0.05882353, 0.04761905, 0.03225806, 0.04347826],
        [0.05882353, 0.        , 0.03333333, 0.02941176, 0.04761905],
        [0.04761905, 0.03333333, 0.        , 0.03571429, 0.02564103],
        [0.03225806, 0.02941176, 0.03571429, 0.        , 0.02325581],
        [0.04347826, 0.04761905, 0.02564103, 0.02325581, 0.        ]])
    x = wpgma(M)
    assert((((0, 1), 4), (2, 3)) == x)


