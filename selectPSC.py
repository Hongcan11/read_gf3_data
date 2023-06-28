import h5py
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi

def selectPSc(stackHDF, D_A_thr=0.25):
    with h5py.File(stackHDF, "a") as f:
        dataset = f['SLC/20230217']
        AOIdataset = dataset[6000:6500, 6000:6500]
    azSize, rgSize = AOIdataset.shape
    beta0array = np.zeros((azSize, rgSize))
    DN = np.abs(AOIdataset)
    Ks = 0.00000000025
    beta0 = (Ks * np.power(np.abs(DN), 2))[0, 0]
    print("beta0: ", beta0)
    beta0array = np.power(np.abs(AOIdataset), 2) / (
            beta0 ** 2)
    print("beta0array ", beta0array, beta0array.shape)
    meanBeta_dB = 10 * np.log10(beta0array + 1e-14)
    print("meanBetaDB: ", meanBeta_dB)
    meanAmp = np.sqrt(np.power(10, (meanBeta_dB / 10)))
    stdAmp = np.nanstd(np.sqrt(beta0array))
    D_A = stdAmp / meanAmp
    print("meanAmp: ", meanAmp)
    print("stdAmp: ", stdAmp)
    print("D_A: ", D_A)
    peaks = peak_local_max(meanAmp)
    print("peaks: ", peaks, peaks.shape)
    azIdx = peaks.T[0]
    rgIdx = peaks.T[1]
    print("azIdx, rgIdx: ", azIdx, rgIdx)
    # use D_A threshold:
    peak_D_A = D_A[peaks.T[0], peaks.T[1]]
    print("peak_D_A: ", peak_D_A, peak_D_A.shape)
    azIdx = azIdx[peak_D_A < D_A_thr]
    rgIdx = rgIdx[peak_D_A < D_A_thr]
    print("select azIdx, rgIdx: ", azIdx, rgIdx)
    print("azidx, rgidx size: ", azIdx.size, rgIdx.size)
    # # read SLC and IFG:
    SLCarray = np.zeros((azIdx.size), dtype="cfloat")
    SLC = AOIdataset[:]
    SLCarray[:] = SLC[azIdx, rgIdx]
    
    # prepare psc dict:
    psc = {
        "dataset": AOIdataset,
        "azIdx": azIdx,
        "rgIdx": rgIdx,
        "SLC": SLCarray,
        "D_A": D_A[azIdx, rgIdx],
    }
    print("D_A: ", psc["D_A"].shape)
    return psc

def peak_local_max(image, min_distance=1, threshold_abs=None,
                   threshold_rel=None, exclude_border=True, indices=True,
                   num_peaks=np.inf, footprint=None, labels=None,
                   num_peaks_per_label=np.inf):

    out = np.zeros_like(image, dtype=np.bool)

    threshold_abs = threshold_abs if threshold_abs is not None else image.min()

    if isinstance(exclude_border, bool):
        exclude_border = (min_distance if exclude_border else 0,) * image.ndim
    elif isinstance(exclude_border, int):
        if exclude_border < 0:
            raise ValueError("`exclude_border` cannot be a negative value")
        exclude_border = (exclude_border,) * image.ndim
    elif isinstance(exclude_border, tuple):
        if len(exclude_border) != image.ndim:
            raise ValueError(
                "`exclude_border` should have the same length as the "
                "dimensionality of the image.")
        for exclude in exclude_border:
            if not isinstance(exclude, int):
                raise ValueError(
                    "`exclude_border`, when expressed as a tuple, must only "
                    "contain ints."
                )
            if exclude < 0:
                raise ValueError(
                    "`exclude_border` cannot contain a negative value")
    else:
        raise TypeError(
            "`exclude_border` must be bool, int, or tuple with the same "
            "length as the dimensionality of the image.")

    # no peak for a trivial image
    if np.all(image == image.flat[0]):
        if indices is True:
            return np.empty((0, image.ndim), np.int)
        else:
            return out

    # In the case of labels, call ndi on each label
    if labels is not None:
        label_values = np.unique(labels)
        # Reorder label values to have consecutive integers (no gaps)
        if np.any(np.diff(label_values) != 1):
            mask = labels >= 1
            labels[mask] = 1 + rank_order(labels[mask])[0].astype(labels.dtype)
        labels = labels.astype(np.int32)

        # create a mask for the non-exclude region
        inner_mask = _exclude_border(np.ones_like(labels, dtype=bool),
                                     exclude_border)

        # For each label, extract a smaller image enclosing the object of
        # interest, identify num_peaks_per_label peaks and mark them in
        # variable out.
        for label_idx, obj in enumerate(ndi.find_objects(labels)):
            img_object = image[obj] * (labels[obj] == label_idx + 1)
            mask = _get_peak_mask(img_object, min_distance, footprint,
                                  threshold_abs, threshold_rel)
            if exclude_border:
                # remove peaks fall in the exclude region
                mask &= inner_mask[obj]
            coordinates = _get_high_intensity_peaks(img_object, mask,
                                                    num_peaks_per_label)
            nd_indices = tuple(coordinates.T)
            mask.fill(False)
            mask[nd_indices] = True
            out[obj] += mask

        if not indices and np.isinf(num_peaks):
            return out

        coordinates = _get_high_intensity_peaks(image, out, num_peaks)
        if indices:
            return coordinates
        else:
            out.fill(False)
            nd_indices = tuple(coordinates.T)
            out[nd_indices] = True
            return out

    # Non maximum filter
    mask = _get_peak_mask(image, min_distance, footprint, threshold_abs,
                          threshold_rel)

    mask = _exclude_border(mask, exclude_border)

    # Select highest intensities (num_peaks)
    coordinates = _get_high_intensity_peaks(image, mask, num_peaks)

    if indices is True:
        return coordinates
    else:
        nd_indices = tuple(coordinates.T)
        out[nd_indices] = True
        return out
    
def rank_order(image):
    flat_image = image.ravel()
    sort_order = flat_image.argsort().astype(np.uint32)
    flat_image = flat_image[sort_order]
    sort_rank = np.zeros_like(sort_order)
    is_different = flat_image[:-1] != flat_image[1:]
    np.cumsum(is_different, out=sort_rank[1:])
    original_values = np.zeros((sort_rank[-1] + 1,), image.dtype)
    original_values[0] = flat_image[0]
    original_values[1:] = flat_image[1:][is_different]
    int_image = np.zeros_like(sort_order)
    int_image[sort_order] = sort_rank
    return (int_image.reshape(image.shape), original_values)

def _exclude_border(mask, exclude_border):
    for i, excluded in enumerate(exclude_border):
        if excluded == 0:
            continue
        mask[(slice(None),) * i + (slice(None, excluded),)] = False
        mask[(slice(None),) * i + (slice(-excluded, None),)] = False
    return mask

def _get_peak_mask(image, min_distance, footprint, threshold_abs,
                   threshold_rel):
    """
    Return the mask containing all peak candidates above thresholds.
    """
    if footprint is not None:
        image_max = ndi.maximum_filter(image, footprint=footprint,
                                       mode='constant')
    else:
        size = 2 * min_distance + 1
        image_max = ndi.maximum_filter(image, size=size, mode='constant')
    mask = image == image_max
    if threshold_rel is not None:
        threshold = max(threshold_abs, threshold_rel * image.max())
    else:
        threshold = threshold_abs
    mask &= image > threshold
    return mask

def _get_high_intensity_peaks(image, mask, num_peaks):
    """
    Return the highest intensity peak coordinates.
    """
    # get coordinates of peaks
    coord = np.nonzero(mask)
    intensities = image[coord]
    # Highest peak first
    idx_maxsort = np.argsort(-intensities)
    coord = np.transpose(coord)[idx_maxsort]
    # select num_peaks peaks
    if len(coord) > num_peaks:
        coord = coord[:num_peaks]
    return coord

def plot_psc(psc):
    D_A = psc["D_A"][:]
    az = psc["azIdx"][:]
    rg = psc["rgIdx"][:]
    print(D_A.shape, D_A)
    fig, ax = plt.subplots(figsize=(10, 10), dpi=200)
    q = ax.scatter(rg, az, c=D_A, s=10)
    ax.set_xlabel("pixels")
    ax.yaxis.set_tick_params(rotation=90)
    ax.set_ylabel("lines")
    ax.set_title("D_A of PSC")
    plt.set_cmap("viridis")
    fig.colorbar(q, ax=ax)
    plt.savefig('/data/tests/hongcan/GF3/inner_mongolia', bbox_inches="tight")

if __name__ == "__main__":
    psc = selectPSc(stackHDF= '/data/tests/hongcan/GF3/cn_inner_mongolia_gf3c_dsc_fsii/gf3.hdf5', D_A_thr=0.25)
    plot_psc(psc)