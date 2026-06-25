import numpy as np

def  normalised_cross_correlation(g1,g2):

    a = (g1 - np.mean(g1)) / (np.std(g1) * len(g1))
    b = (g2 - np.mean(g2)) / (np.std(g2))
    
    cc = np.correlate(a, b, 'full')
    max_cc = np.max(cc)
    min_cc = np.min(cc)

    return max_cc,min_cc