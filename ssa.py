import numpy as np

def running_mean(x,window):
    cumsum = np.cumsum(np.insert(x,0,0))
    return (cumsum[window:]-cumsum[:-window])/window

def SSA(Y,L,period_groups):
    T = Y.size
    assert L <= T/2
    K = T - L + 1

    # Form the trajectory matrix and find the eigen decomp
    X = np.zeros((L,K))
    for i in range(K): X[:,i] = Y[i:(i+L)]
    lamda,P = np.linalg.eig(np.dot(X,X.T))

    # Find the dominant frequency of each eigenvector
    f   = np.zeros(lamda.size)
    fs  = np.fft.fftfreq(f.size,1.)
    ix  = np.argsort(fs)
    fs  = fs[ix]
    eps = 0.99*(fs[1]-fs[0])
    for i in range(f.size):
        ps = np.abs(np.fft.fft(P[:,i]))**2
        ps = ps[ix]
        f[i] = fs[ps.argmax()]
    f = np.abs(f)

    # convert periodicity into frequency
    fgroups = 1/np.asarray(period_groups,dtype=float)
    fgroups = np.hstack([0,fgroups])
    
    # Build an approximation of X by taking a subset of the
    # decomposition. This approximation is formed by taking
    # eigenvectors whose dominant frequency is close to the targetted
    # values.
    Xt = np.zeros((fgroups.size,)+X.shape)
    P = P.real #bug fix for sm 
    for i in range(f.size):
        g = np.where(np.abs(fgroups-f[i]) < eps)[0]
        if g.size == 0: continue
        Xt[g[0]] += np.dot(np.outer(P[:,i],P[:,i]),X)

    # Now we reconstruct the signal by taking a mean of all the
    # approximations.
    Yt = np.zeros((fgroups.size,Y.size))
    c  = np.zeros((fgroups.size,Y.size))
    for g in range(fgroups.size):
        for i in range(K): 
            Yt[g,i:(i+L)] += Xt[g,:,i]
            c [g,i:(i+L)] += 1
    Yt /= c
    
    return Yt

def GetResidual(Y,window=120,groups=[12,6,4,3,2]):
    """Extracts anomalies from a signal.

    Parameters
    ----------
    Y : numpy.ndarray
        the one-dimensional, non-masked array of data
    window : int
        the window size, the maximum periodicity you wish to search
        for in the input signal
    groups: array-like of ints
        a list of data periodicities to extract from the input signal

    Returns
    -------
    R : numpy.ndarray
        the signal minus the trend and the group periodicities
        specified

    """
    Yt = SSA(Y,window,groups)
    R  = Y-Yt.sum(axis=0)
    return R

def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index."""

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero() 

    # We need to start things after the change in "condition". Therefore, 
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size] # Edit

    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx
    
if __name__ == "__main__":
    from netCDF4 import Dataset
    import pylab as plt

    dset = Dataset("/work/ILAMB/MODELS/CLM/CLM50r243GSWP3/gpp_Lmon_CLM50r243GSWP3_historical_r1i1p1_185001-201012.nc")
    i    = np.abs(dset.variables["lat"][...]+  5.).argmin()
    j    = np.abs(dset.variables["lon"][...]-300.).argmin()
    t    = dset.variables["time"][...  ]/365.+1850.
    gpp  = dset.variables["gpp" ][:,i,j]
    groups = [12,6,4,3,2]
    decomp = SSA(gpp,120,groups)
    
    fig,axs = plt.subplots(nrows=len(groups)+2,tight_layout=True)
    axs[0].plot(t,gpp,'-')
    axs[0].set_ylim(gpp.min(),gpp.max())
    axs[1].set_ylim(gpp.min(),gpp.max())
    for g in range(len(groups)):
        axs[g+1].plot(t,decomp[g],'-')
    axs[-1].plot(t,gpp-decomp.sum(axis=0),'-')
    axs[-1].plot(t,GetResidual(gpp),'--k')
    plt.show()
