from typing import Literal, Tuple

import numpy as np
import numpy.typing as npt


def enrichment_factor_score(
    y_true: npt.ArrayLike,
    y_score: npt.ArrayLike,
    fraction: float = 0.01,
    tie_handling: Literal["strict", "expand"] = "strict",
) -> Tuple[float, float]:
    """
    Calculate the Enrichment Factor (EF), with an option to handle tied scores
    and return the actual fraction of samples used if ties are expanded.

    Parameters
    ----------
    y_true : npt.ArrayLike of shape (n_samples,)
        Ground truth labels for each sample (1 = active, 0 = inactive).
    y_score : npt.ArrayLike of shape (n_samples,)
        Model scores or predicted scores (higher = more likely to be active).
    fraction : float, default=0.01
        The fraction (0 < fraction <= 1) of the total samples to consider as "top."
        For example, fraction=0.01 corresponds to the top 1%.
    tie_handling : {"strict", "expand"}, default="strict"
        "strict" : Select exactly floor(fraction * n_samples) samples.
        "expand" : If the sample at the boundary ties with the next sample's score,
                   include all samples with that same score. This can increase the
                   effective fraction above the user-specified fraction.

    Returns
    -------
    ef_value : float
        The Enrichment Factor (EF) at the specified fraction (or expanded if ties).
    actual_fraction : float
        The actual fraction of samples used.

    Raises
    ------
    ValueError
        - If fraction is not in (0, 1].
        - If fraction is too small such that the top count is 0.
        - If y_true and y_score have different lengths.
    """

    # Convert inputs to numpy arrays
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

    # Sort indices by score in descending order
    sort_indices = np.argsort(-y_score_array)
    # For "strict", we just take the top_n
    if tie_handling == "strict":
        final_top_indices = sort_indices[:top_n]

    # For "expand", we include all samples that tie with the last selected score
    else:  # tie_handling == "expand"
        # First, select top_n
        temp_top_indices = sort_indices[:top_n]
        # Score of the boundary sample
        boundary_score = y_score_array[temp_top_indices[-1]]
        # Check if there are additional samples (after top_n) with the same score
        tie_indices = []
        for idx in sort_indices[top_n:]:
            if y_score_array[idx] == boundary_score:
                tie_indices.append(idx)
            else:
                break  # because we are in descending order, once we pass the boundary score, we can stop

        final_top_indices = np.concatenate([temp_top_indices, tie_indices]) if tie_indices else temp_top_indices

    # Count how many active samples are in the final top range
    top_actives = np.sum(y_true_array[final_top_indices])
    final_top_count = len(final_top_indices)  # may exceed top_n if we expanded

    # Compute Enrichment Factor = (active rate in top fraction) / (active rate overall)
    top_active_ratio = top_actives / final_top_count
    overall_active_ratio = total_actives / n_samples if n_samples else 0

    # If there are no active samples overall, EF cannot be defined. Return 0.0 or handle as appropriate.
    if overall_active_ratio == 0:
        ef_value = 0.0
    else:
        ef_value = top_active_ratio / overall_active_ratio

    # Compute the actual fraction used (might exceed 'fraction' in "expand" mode)
    actual_fraction = final_top_count / n_samples

    return ef_value, actual_fraction


if __name__ == "__main__":
    # Example data with some tied scores
    y_true_example = [1, 1, 0, 1, 0, 1, 0, 0]
    y_score_example = [0.80, 0.80, 0.78, 0.78, 0.55, 0.55, 0.30, 0.10]

    print("Scores:", y_score_example)

    # Strict approach: top 25% => top 2 samples
    ef_strict, frac_strict = enrichment_factor_score(
        y_true_example, y_score_example, fraction=0.25, tie_handling="strict"
    )
    print("[strict] EF(25%):", ef_strict)
    print("[strict] Actual fraction used:", frac_strict)

    # Expand approach: top 25% => might include more if there's a tie at the boundary
    ef_expand, frac_expand = enrichment_factor_score(
        y_true_example, y_score_example, fraction=0.25, tie_handling="expand"
    )
    print("[expand] EF(25%):", ef_expand)
    print("[expand] Actual fraction used:", frac_expand)
