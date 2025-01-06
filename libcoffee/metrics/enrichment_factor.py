import numpy as np
import numpy.typing as npt


def enrichment_factor_score(y_true: npt.ArrayLike, y_score: npt.ArrayLike, fraction: float = 0.01) -> float:
    """
    Calculate the Enrichment Factor (EF).

    Parameters
    ----------
    y_true : npt.ArrayLike of shape (n_samples,)
        Ground truth labels for each sample, where 1 indicates "active" and 0 indicates "inactive".
    y_score : npt.ArrayLike of shape (n_samples,)
        Model scores or predicted scores for each sample. A higher score indicates a higher likelihood of being active.
    fraction : float, default=0.01
        The fraction (between 0 and 1) of the total samples considered as the "top" range.
        For instance, fraction=0.01 corresponds to the top 1%.

    Returns
    -------
    float
        The Enrichment Factor (EF) at the specified fraction.

    Raises
    ------
    ValueError
        - If fraction is not in (0, 1].
        - If fraction is too small such that the top count is 0.
        - If y_true and y_score have different lengths.
    """

    # Convert inputs to NumPy arrays
    y_true_array = np.asarray(y_true)
    y_score_array = np.asarray(y_score)

    # Input validation
    if not (0 < fraction <= 1):
        raise ValueError("fraction must be between 0 and 1.")

    if y_true_array.shape[0] != y_score_array.shape[0]:
        raise ValueError("y_true and y_score must have the same length.")

    # Total number of samples and total number of active samples
    n_samples = y_true_array.shape[0]
    total_actives = np.sum(y_true_array)

    # Determine the number of top samples based on fraction
    top_n = int(n_samples * fraction)
    if top_n < 1:
        raise ValueError("The fraction is too small; no top samples are selected. Please increase the fraction.")

    # Sort indices by score in descending order, then select the top_n
    sort_indices = np.argsort(-y_score_array)  # negative sign to sort in descending order
    top_indices = sort_indices[:top_n]  # TODO: 同率の場合の処理を追加する

    # Count how many active samples are in the top range
    top_actives = np.sum(y_true_array[top_indices])

    # Compute Enrichment Factor = (active rate in top fraction) / (active rate overall)
    top_active_ratio = top_actives / top_n
    overall_active_ratio = total_actives / n_samples

    # If there are no active samples overall, EF cannot be defined. Return 0.0 or handle as appropriate.
    if overall_active_ratio == 0:
        return 0.0

    ef_value = top_active_ratio / overall_active_ratio
    return float(ef_value)


# Example usage:
if __name__ == "__main__":
    # Example data
    y_true_example = [1, 1, 0, 1, 0, 1, 0, 0]
    y_score_example = [0.95, 0.80, 0.78, 0.60, 0.55, 0.40, 0.30, 0.10]

    # Compute EF at top 10%
    ef_10 = enrichment_factor_score(y_true_example, y_score_example, fraction=0.1)
    print("EF(10%):", ef_10)

    # Compute EF at top 25%
    ef_25 = enrichment_factor_score(y_true_example, y_score_example, fraction=0.25)
    print("EF(25%):", ef_25)
