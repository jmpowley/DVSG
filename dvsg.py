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
        The normalised velocity map, with values in the range [-1, 1].
        If all finite values in velocity_map are identical, returns NaNs.
    """

    min_val, max_val = np.nanmin(velocity_map), np.nanmax(velocity_map)

    if min_val == max_val:
        return np.full_like(velocity_map, np.nan)  # Avoid division by zero

    normalised_velocity_map = 2 * (velocity_map - min_val) / (max_val - min_val) - 1

    return normalised_velocity_map

def denormalise_velocity_map(normalised_velocity_map: np.ndarray, max_velocity, min_velocity):
    """Denormalises a velocity map from the range [-1, 1] back to its original scale.

    This function requires the minimum and maximum of the original velocity map before normalisation.

    The denormalisation formula used is:

        x = ((x' + 1) / 2) * (max_velocity - min_velocity) + min_velocity

    Parameters
    ----------
    normalised_velocity_map : np.ndarray
        The normalised velocity map in the range [-1, 1].
    max_velocity : float/int
        The maximum value of the original unnormalised velocity map.
    min_velocity : float/int
        The minimum value of the original unnormalised velocity map.

    Returns
    -------
    unnormalised_velocity_map : np.ndarray
        The unnormalised velocity map, with values in the range [min_velocity, max_velocity].
        If all finite values in normalised_velocity_map are identical, returns NaNs.
    """

    if not type(max_velocity) is int or float:
        raise Exception("max_velocity must be an int or float")  # Ensure max_velocity is float/int
    if not type(min_velocity) is int or float:
        raise Exception("max_velocity must be an int or float")  # Ensure min_velocity is float/int
    if min_velocity == max_velocity:
        return np.full_like(normalised_velocity_map, np.nan)  # Avoid division by zero

    unnormalised_velocity_map = ((normalised_velocity_map + 1) / 2) * (max_velocity - min_velocity) + min_velocity

    return unnormalised_velocity_map

def calculate_DVSG(sv_map: np.ndarray, gv_map: np.ndarray):
    """Calculates the DVSG value for a galaxy's stellar and gas velocity map. 
    These maps are assumed to be normalised between -1 and 1.

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
        The DVSG value for the galaxy (sum of all values in dvsg_map, normalised by
        the number of finite elements).
    """

    if not np.shape(sv_map) == np.shape(gv_map):
        raise Exception("sv_map and gv_map must have the same shape. Currently have shapes " + str(np.shape(sv_map)) + " and " + str(np.shape(gv_map)))
    
    dvsg_map = np.abs(sv_map - gv_map)
    dvsg = np.nansum(dvsg_map) / np.count_nonzero(np.isfinite(dvsg_map))
    
    return dvsg_map, dvsg