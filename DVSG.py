import numpy as np

def exclude_above_five_sigma(velocity_map: np.ndarray):
    """Excludes any values in a velocity map greater than 5 standard deviations
    from the mean velocity of the map.
    
    NaN values are ignored when computing the mean and standard deviation.

    The standard deviation used is the sample standard deviation (ddof=1).

    Parameters
    ----------
    velocity_map : np.ndarray
        The original velocity map.

    Returns
    -------
    excluded_velocity_map : np.ndarray
        The original velocity map, but with all velocity values greater 
        than 5 standard deviations from the mean set to NaN.
    """

    velocity_above_five_sigma = np.nanmean(velocity_map) + 5 * np.nanstd(velocity_map, ddof=1)
    velocity_below_five_sigma = np.nanmean(velocity_map) - 5 * np.nanstd(velocity_map, ddof=1)
    
    excluded_velocity_map = velocity_map.copy()

    excluded_velocity_map[velocity_map > velocity_above_five_sigma] = np.nan
    excluded_velocity_map[velocity_map < velocity_below_five_sigma] = np.nan

    return excluded_velocity_map

def normalise_velocity_map(velocity_map: np.ndarray):
    """Normalises the given velocity map to the range [-1, 1].

    Normalisation uses the formula:

        x' = 2 * ((x - min(x)) / (max(x) - min(x))) - 1 

    NaN values are ignored when computing the minimum and maximum.

    Parameters
    ----------
    velocity_map : np.ndarray
        The unnormalised velocity map.

    Returns
    -------
    normalised_velocity_map : np.ndarray
        The normalized velocity map, with values in the range [-1, 1].
        If all finite values in velocity_map are identical, returns NaNs.
    """

    min_val, max_val = np.nanmin(velocity_map), np.nanmax(velocity_map)

    if min_val == max_val:
        return np.full_like(velocity_map, np.nan)  # Avoid division by zero

    normalised_velocity_map = 2 * (velocity_map - min_val) / (max_val - min_val) - 1

    return normalised_velocity_map

def denormalise_velocity_map(normalized_velocity, original_velocity):
    # Denormalize values back to original range
    norm_min = -1
    norm_max = 1
    return (normalized_velocity - norm_min) / (norm_max - norm_min) * (np.nanmax(original_velocity) - np.nanmin(original_velocity)) + np.nanmin(original_velocity)

def calculate_DVSG(sv_map, gv_map):
    """Calculates the DVSG value for a galaxy's stellar and gas velocity map. 
    These maps are assumed to be normalized between -1 and 1.

    Parameters
    ----------
    sv_map : np.ndarray
        The stellar velocity map.
    gv_map : np.ndarray
        The gas velocity map.

    sv_map and gv_map must have the same shape.

    Returns
    -------
    dvsg_map : np.ndarray
        The absolute difference in the stellar and gas velocity maps.
    dvsg : float
        The DVSG value for the galaxy (sum of all values in dvsg_map, normalized by the number of finite elements).
    """
    
    dvsg_map = np.abs(sv_map - gv_map)
    dvsg = np.nansum(dvsg_map) / np.count_nonzero(np.isfinite(dvsg_map))
    
    return dvsg_map, dvsg