import pandas as pd
import random
from rapidfuzz import process, fuzz
import time
import numpy as np
from scipy.stats import mode
from sklearn.metrics import pairwise_distances
import spoa
from numba import njit, prange
from rapidfuzz import fuzz
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram


# --- 2. DATA GENERATION & UTILS ---
def mutate_sequence(seq, error_rate=0.10):
    if error_rate <= 0:
        return seq

    seq_list = list(seq)
    new_seq = []
    for base in seq_list:
        if random.random() < error_rate:
            r = random.random()
            if r < 0.5:
                new_seq.append(random.choice("ACGT"))
            elif r < 0.75:
                new_seq.append(base)
                new_seq.append(random.choice("ACGT"))
            else:
                pass
        else:
            new_seq.append(base)
    return "".join(new_seq)


def gen_dna(k): return "".join(random.choices("ACGT", k=k))


# (Generating fresh data to ensure variables exist for the loop below)
n_items = 3_000
pool_bc = [gen_dna(12) for _ in range(n_items)]
pool_ins = [gen_dna(random.randint(300, 500)) for _ in range(n_items)]
n_repeat_bc = 500
pool_bc += pool_bc[0:n_repeat_bc]
pool_ins += [gen_dna(random.randint(300, 500)) for _ in range(n_repeat_bc)]
pool_bc += pool_bc[0:n_repeat_bc]
pool_ins += [gen_dna(random.randint(300, 500)) for _ in range(n_repeat_bc)]

# Add very similar sequence that will defeat clustering!
pool_bc += [mutate_sequence(pool_bc[i], 0.03) for i in range(n_repeat_bc)]
pool_ins += [mutate_sequence(pool_ins[i], 0.02) for i in range(n_repeat_bc)]

data = []
for i in range(len(pool_bc)):
    n_reads = random.randint(5, 20)
    for _ in range(n_reads):
        data.append({
            "ID": i,
            "Barcode": mutate_sequence(pool_bc[i], 0.01),
            "Insert": mutate_sequence(pool_ins[i], 0.01)
        })

df = pd.DataFrame(data)


# Assuming rapidfuzz or fuzzywuzzy is available based on context


# ==========================================
# TIMING UTILITY
# ==========================================
class GlobalTimer:
    def __init__(self):
        self.timings = {}
        self.starts = {}

    def start(self, name):
        if name not in self.starts:  # Don't overwrite if nested recursion causes restart
            self.starts[name] = time.perf_counter()

    def stop(self, name):
        if name in self.starts:
            elapsed = time.perf_counter() - self.starts[name]
            self.timings[name] = self.timings.get(name, 0) + elapsed
            del self.starts[name]

    def report(self):
        print("\n" + "=" * 40)
        print(f"{'OPERATION':<30} | {'TIME (s)':<10}")
        print("-" * 43)
        for name, val in sorted(self.timings.items(), key=lambda x: x[1], reverse=True):
            print(f"{name:<30} | {val:.4f}")
        print("=" * 40 + "\n")


timer = GlobalTimer()


# ==========================================
# CORE FUNCTIONS
# ==========================================
def remove_invariant_columns(msa_int):
    """
    Identifies and removes columns where values do not change (invariant).
    Returns the reduced matrix and the indices of the kept columns.
    """
    if msa_int.size == 0:
        return msa_int, np.array([])

    mins = msa_int.min(axis=0)
    maxs = msa_int.max(axis=0)

    # If min == max, the column has only one unique value (invariant)
    is_variable = (mins != maxs)
    kept_indices = np.where(is_variable)[0]

    msa_reduced = msa_int[:, kept_indices]
    return msa_reduced, kept_indices


def mutate_integers_fast(seed_int_arr, error_rate, n_copies, mutation_choices):
    """
    Mutates an integer array directly using NumPy.
    Orders of magnitude faster than string manipulation.
    """
    L = len(seed_int_arr)

    # 1. Tile the seed (Create N copies)
    copies = np.tile(seed_int_arr, (n_copies, 1))

    # 2. Generate Error Mask
    mask = np.random.random((n_copies, L)) < error_rate

    # 3. Apply Mutations
    n_mutations = np.sum(mask)
    if n_mutations > 0:
        replacements = np.random.choice(mutation_choices, size=n_mutations)
        copies[mask] = replacements

    return copies


def encode_msa(msa_strings, alphabet="ACGTN-"):
    timer.start("encode_msa")
    """
    Vectorized encoding using ASCII byte views.
    Assumes inputs are ASCII (standard for DNA).
    """
    N = len(msa_strings)
    if N == 0:
        timer.stop("encode_msa")
        return np.array([]), len(alphabet)

    max_len = max(len(s) for s in msa_strings)

    # 1. Pad and convert to single byte-string buffer
    padded_block = "".join([s.ljust(max_len, '-') for s in msa_strings]).encode('ascii')

    # 2. View as NumPy int8 array directly
    arr_view = np.frombuffer(padded_block, dtype=np.int8).copy()
    arr_view = arr_view.reshape(N, max_len)

    # 3. Fast Translation Table
    lookup = np.zeros(256, dtype=np.int8) + (len(alphabet) - 1)
    for idx, char in enumerate(alphabet):
        lookup[ord(char)] = idx

    # 4. Apply translation
    msa_int = lookup[arr_view]

    timer.stop("encode_msa")
    return msa_int, len(alphabet)


def get_marginals(msa_int, vocab_size):
    # Create on hot matrix
    # Use float32 for faster matrix multiplication later
    one_hot = np.eye(vocab_size, dtype=np.float32)[msa_int]
    # Calculate frequencies of each base/N/- at each position in MSA
    P_i = one_hot.mean(axis=0)
    return one_hot, P_i


@njit(parallel=True, fastmath=True)
def compute_mi_from_joint(joint_probs_flat, P_i_flat, L, A):
    """
    Computes MI from the pre-calculated joint probability matrix.
    joint_probs_flat: (L*A, L*A) matrix containing P(x,y)
    P_i_flat: (L*A) vector containing P(x)
    """
    mi = np.zeros((L, L), dtype=np.float32)
    epsilon = 1e-12
    # Parallel loop over positions
    for i in prange(L):
        for j in range(i + 1, L):  # Triangle only (symmetric)
            mi_val = 0.0

            # Loop over alphabet characters
            for a in range(A):
                p_x = P_i_flat[i * A + a]

                for b in range(A):
                    p_y = P_i_flat[j * A + b]

                    # Look up joint prob directly from flat matrix
                    row_idx = i * A + a
                    col_idx = j * A + b
                    p_xy = joint_probs_flat[row_idx, col_idx]

                    if p_xy > epsilon:
                        product = p_x * p_y
                        if product > epsilon:
                            mi_val += p_xy * np.log(p_xy / product)

            mi[i, j] = mi_val
            mi[j, i] = mi_val
    return mi


def get_elbow_columns(mi_matrix, exclusion_distance=5, verbose=0):
    timer.start("get_elbow_columns")
    L = mi_matrix.shape[0]
    # Note: When using reduced matrix, exclusion_distance refers to indices in the reduced set
    mask = np.triu(np.ones((L, L), dtype=bool), k=exclusion_distance + 1)
    rows, cols = np.where(mask)
    scores = mi_matrix[rows, cols]

    if len(scores) == 0:
        timer.stop("get_elbow_columns")
        return np.array([])
    # 1. Sort Descending
    sorted_indices = np.argsort(scores)[::-1]
    sorted_scores = scores[sorted_indices]

    # 2. Determine "Signal End" (Noise Floor Truncation)
    noise_floor = np.percentile(scores, 25)
    valid_mask = sorted_scores > noise_floor
    n_signal_points = np.sum(valid_mask)
    min_points = min(len(sorted_scores), 5)
    cutoff_idx = max(n_signal_points, min_points)

    # 3. Truncate for Geometry Calculation
    curve_y = sorted_scores[:cutoff_idx]
    n_points = len(curve_y)

    if n_points < 3:
        timer.stop("get_elbow_columns")
        return np.unique(np.concatenate([rows[sorted_indices[:n_points]], cols[sorted_indices[:n_points]]]))
    # 4. Standard Kneedle on Truncated Curve
    x_norm = np.linspace(0, 1, n_points)
    y_norm = (curve_y - curve_y.min()) / (curve_y.max() - curve_y.min() + 1e-9)

    line_vec = np.array([1.0, -1.0])
    line_vec = line_vec / np.linalg.norm(line_vec)
    vec_from_start = np.stack([x_norm, y_norm - 1.0], axis=1)
    distances = np.abs(vec_from_start[:, 0] * line_vec[1] - vec_from_start[:, 1] * line_vec[0])

    # 5. Select
    elbow_idx = np.argmax(distances)
    if elbow_idx == 0:
        for i in range(1, n_points - 1):
            if curve_y[i] >= (curve_y[0] * 0.95):
                elbow_idx = i
            else:
                break
    n_selected = elbow_idx + 1
    selected_indices = sorted_indices[:n_selected]
    selected_columns = np.unique(np.concatenate([rows[selected_indices], cols[selected_indices]]))

    timer.stop("get_elbow_columns")
    return selected_columns


def plot_dendrogram(model, **kwargs):
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count
    linkage_matrix = np.column_stack([model.children_, model.distances_, counts]).astype(float)
    dendrogram(linkage_matrix, **kwargs)


def fast_consensus(msa_list, min_threshold=0.5, high_conf_threshold=0.85):
    """
    Vectorized consensus with strict ambiguity detection.

    Rules:
      1. TIE (e.g., 2 A's, 2 G's) -> 'N'
      2. Support < min_threshold  -> 'N'
      3. Support < high_conf      -> Lower Case (e.g., 'a')
      4. Support >= high_conf     -> Upper Case (e.g., 'A')
    """
    if not msa_list:
        return ""

    # 1. Setup Vectorized Arrays
    arr = np.array([list(s) for s in msa_list], dtype='S1')
    arr_int = arr.view(np.uint8)
    n_seqs, length = arr_int.shape

    GAP = 45  # '-'
    N_VAL = 78  # 'N'
    TO_LOWER = 32  # ASCII offset for A->a

    # 2. Count Frequencies (Vectorized)
    unique_chars = np.unique(arr_int)
    # Shape: (num_unique, length)
    counts = (arr_int == unique_chars[:, None, None]).sum(axis=1)

    # 3. Determine Winners & Detect Ties
    max_counts = counts.max(axis=0)
    max_indices = counts.argmax(axis=0)

    # Check for ties: Count how many rows in 'counts' equal the 'max_count'
    # If sum > 1, it means two or more bases are tied for the winner.
    tie_counts = (counts == max_counts).sum(axis=0)
    is_tie = tie_counts > 1

    # Map to Consensus Candidates
    consensus_arr = unique_chars[max_indices]
    freqs = max_counts / n_seqs

    # 4. Apply Rules
    is_gap = (consensus_arr == GAP)

    # Rule A: Ties or Low Support -> 'N'
    # Note: We enforce N if it's a tie, OR if freq is below min threshold.
    ambiguous_mask = (is_tie) | (freqs < min_threshold)

    # Only overwrite non-gaps (preserving gaps allows us to strip them later)
    # If a gap is the clear winner, we leave it as gap.
    # If a gap is TIED with a base (e.g. 2 gaps, 2 As), is_tie=True, so it becomes N.
    consensus_arr[ambiguous_mask & (~is_gap)] = N_VAL

    # Rule B: Mid Confidence -> Lower Case
    # Valid only if NOT ambiguous, NOT a gap, and freq < high_threshold
    mid_conf_mask = (~ambiguous_mask) & (freqs < high_conf_threshold) & (~is_gap)

    # Safety: Only lowercase actual letters (A-Z)
    is_upper_alpha = (consensus_arr >= 65) & (consensus_arr <= 90)
    consensus_arr[mid_conf_mask & is_upper_alpha] += TO_LOWER

    # 5. Output
    # Remove gaps
    final_int_seq = consensus_arr[consensus_arr != GAP]
    return final_int_seq.tobytes().decode('utf-8')


def get_poa_consensus(barcode_series, bias_towards=None):
    timer.start("SPOA_Consensus_Total")
    barcode_list = barcode_series.tolist()
    if bias_towards is not None:
        barcode_list.append(bias_towards)

    _, msa = spoa.poa(barcode_list, algorithm=1)
    consensus = fast_consensus(msa)
    timer.stop("SPOA_Consensus_Total")
    return consensus


# ==========================================
# OPTIMIZED MI COMPUTATION
# ==========================================
@njit(parallel=True, fastmath=True)
def compute_mi_sparse(joint_probs_flat, P_i_flat, L, A, position_pairs):
    """
    Computes MI only for specified position pairs.
    Much faster when you only need a subset of interactions.

    position_pairs: Nx2 array of (i, j) pairs to compute
    """
    n_pairs = position_pairs.shape[0]
    mi_values = np.zeros(n_pairs, dtype=np.float32)
    epsilon = 1e-12
    for pair_idx in prange(n_pairs):
        i = position_pairs[pair_idx, 0]
        j = position_pairs[pair_idx, 1]

        mi_val = 0.0
        for a in range(A):
            p_x = P_i_flat[i * A + a]

            for b in range(A):
                p_y = P_i_flat[j * A + b]

                row_idx = i * A + a
                col_idx = j * A + b
                p_xy = joint_probs_flat[row_idx, col_idx]

                if p_xy > epsilon:
                    product = p_x * p_y
                    if product > epsilon:
                        mi_val += p_xy * np.log(p_xy / product)

        mi_values[pair_idx] = mi_val

    return mi_values


@njit(parallel=True, fastmath=True)
def compute_mi_from_joint_optimized(joint_probs_flat, P_i_flat, L, A):
    """
    Optimized version with better memory access patterns.
    """
    mi = np.zeros((L, L), dtype=np.float32)
    epsilon = 1e-12
    # Precompute log(P_i) to avoid redundant calculations
    log_P_i = np.zeros(L * A, dtype=np.float32)
    for idx in range(L * A):
        if P_i_flat[idx] > epsilon:
            log_P_i[idx] = np.log(P_i_flat[idx])
    for i in prange(L):
        for j in range(i + 1, L):
            mi_val = 0.0

            for a in range(A):
                p_x = P_i_flat[i * A + a]
                if p_x <= epsilon:
                    continue

                log_p_x = log_P_i[i * A + a]

                for b in range(A):
                    p_y = P_i_flat[j * A + b]
                    if p_y <= epsilon:
                        continue

                    row_idx = i * A + a
                    col_idx = j * A + b
                    p_xy = joint_probs_flat[row_idx, col_idx]

                    if p_xy > epsilon:
                        log_p_y = log_P_i[j * A + b]
                        mi_val += p_xy * (np.log(p_xy) - log_p_x - log_p_y)

            mi[i, j] = mi_val
            mi[j, i] = mi_val

    return mi


def compute_mi_scores_optimized_v2(one_hot, P_i, timer):
    """
    Enhanced version with precomputed logs.
    """
    timer.start("compute_mi_hybrid_execution")

    N, L, A = one_hot.shape

    # BLAS multiplication (unchanged - already optimal)
    flat_view = one_hot.reshape(N, L * A)
    joint_counts = np.dot(flat_view.T, flat_view)
    joint_probs_flat = joint_counts / N

    # Use optimized Numba function
    result = compute_mi_from_joint_optimized(joint_probs_flat, P_i.flatten(), L, A)

    timer.stop("compute_mi_hybrid_execution")
    return result


# ==========================================
# VECTORIZED MAX ROW RATIO
# ==========================================
def calc_max_row_ratio_vectorized(mi_matrix_reduced, n_total_cols, top_pct=0.01, bottom_pct=0.50):
    """
    Fully vectorized version - removes all Python loops.
    ~10-20x faster for typical sizes.
    """
    R = mi_matrix_reduced.shape[0]
    if R == 0:
        return 0.0

    C = mi_matrix_reduced.shape[1]
    if C == 0:
        return 0.0

    # Sort once for all rows
    sorted_rows = np.sort(mi_matrix_reduced, axis=1)

    n_zeros = n_total_cols - C
    k_top = max(1, int(n_total_cols * top_pct))
    k_bottom = max(1, int(n_total_cols * bottom_pct))

    # --- Vectorized Top Mean ---
    if k_top <= C:
        top_mean = np.mean(sorted_rows[:, -k_top:], axis=1)
    else:
        sum_top = np.sum(sorted_rows, axis=1)
        top_mean = sum_top / k_top

    # --- Vectorized Bottom Mean ---
    if k_bottom <= n_zeros:
        bottom_mean = np.zeros(R, dtype=np.float32)
    else:
        n_needed = k_bottom - n_zeros
        sum_bottom = np.sum(sorted_rows[:, :n_needed], axis=1)
        bottom_mean = sum_bottom / k_bottom

    # Vectorized ratio calculation
    damping = 0.01
    ratios = (top_mean + damping) / (bottom_mean + damping)

    return np.max(ratios)


# ==========================================
# OPTIMIZED INTEGRITY CHECK
# ==========================================


# ==========================================
# ADDITIONAL OPTIMIZATION: BATCH PROCESSING
# ==========================================
@njit(parallel=True)
def mutate_integers_fast_batch(seed_int_arr, error_rate, n_copies_per_batch, n_batches, mutation_choices):
    """
    Generate multiple batches of mutations at once for better cache utilization.
    """
    L = len(seed_int_arr)
    total_copies = n_copies_per_batch * n_batches

    # Preallocate entire output
    all_copies = np.empty((total_copies, L), dtype=seed_int_arr.dtype)

    # Fill in parallel
    for batch_idx in prange(n_batches):
        start_idx = batch_idx * n_copies_per_batch
        end_idx = start_idx + n_copies_per_batch

        # Tile seed
        for i in range(n_copies_per_batch):
            all_copies[start_idx + i] = seed_int_arr

        # Generate mutations for this batch
        batch_size = n_copies_per_batch * L
        mask = np.random.random(batch_size) < error_rate

        if np.sum(mask) > 0:
            replacements = np.random.choice(mutation_choices, size=np.sum(mask))

            # Apply mutations
            flat_view = all_copies[start_idx:end_idx].ravel()
            flat_view[mask] = replacements

    return all_copies


# ==========================================
# 1. PRECOMPUTE MI LOOKUP TABLE
# ==========================================
@njit
def create_mi_lookup(n_members):
    """
    Creates a 3D lookup table for MI terms.
    Dimensions: [count_xy, count_x, count_y]
    Pre-calculates: (N_xy/N) * log((N_xy*N) / (N_x*N_y))
    Access is O(1) integer lookup.
    """
    # Max possible count is n_members.
    # We add +1 for 0-based indexing.
    lut = np.zeros((n_members + 1, n_members + 1, n_members + 1), dtype=np.float32)

    N = float(n_members)

    for c_xy in range(1, n_members + 1):
        for c_x in range(1, n_members + 1):
            for c_y in range(1, n_members + 1):
                # Valid counts only (count_xy cannot exceed x or y)
                if c_xy <= c_x and c_xy <= c_y:
                    p_xy = c_xy / N
                    p_x = c_x / N
                    p_y = c_y / N

                    term = p_xy * np.log(p_xy / (p_x * p_y))
                    lut[c_xy, c_x, c_y] = term

    return lut


# ==========================================
# 2. THE CORE KERNEL (Integer Only)
# ==========================================
@njit(fastmath=True)
def _fast_score_kernel(msa, n_members, L, variant_indices, n_var, lut, top_pct, bottom_pct):
    """
    Computes Max Row Ratio using:
    - Integer counting (No OneHot, No Dot Product)
    - Lookup Table (No Log, No Div)
    """
    # 1. Compute Joint Counts & Marginals
    # We use a flat buffer for the counts.
    # Size: n_var * 5 (ACGTN)
    # We only care about columns mapped in variant_indices

    # vocab size fixed to 5 (ACGTN) for speed
    V = 5

    # Counts: (n_var, V)
    marginals = np.zeros((n_var, V), dtype=np.int32)

    # Joint: (n_var, n_var, V, V)
    # For small n_var (<50), this fits in cache.
    # We only fill upper triangle.
    joint = np.zeros((n_var, n_var, V, V), dtype=np.int32)

    for r in range(n_members):
        # We iterate variants, not full L
        for i in range(n_var):
            col_i = variant_indices[i]
            char_i = msa[r, col_i]
            if char_i >= V: continue  # Skip gaps/invalid if any

            marginals[i, char_i] += 1

            # Cross terms
            for j in range(i + 1, n_var):
                col_j = variant_indices[j]
                char_j = msa[r, col_j]
                if char_j >= V: continue

                joint[i, j, char_i, char_j] += 1

    # 2. Compute Score using LUT
    max_ratio = 0.0
    damping = 0.01

    # Limits
    k_top = int(L * top_pct)
    if k_top < 1: k_top = 1
    k_bottom = int(L * bottom_pct)
    if k_bottom < 1: k_bottom = 1
    n_zeros_implicit = L - n_var

    # Reusable row buffer
    row_mis = np.zeros(n_var, dtype=np.float32)

    for i in range(n_var):
        # Calculate MI row i vs all j
        for j in range(n_var):
            if i == j:
                row_mis[j] = 0.0
                continue

            # Sort indices for Upper Triangle Access
            if i < j:
                u, v = i, j
            else:
                u, v = j, i

            mi_val = 0.0

            # Sum 5x5 alphabet
            for a in range(V):
                c_x = marginals[i, a]
                if c_x == 0: continue

                for b in range(V):
                    c_y = marginals[j, b]
                    if c_y == 0: continue

                    # Access Joint (always u, v)
                    # Mapping: if i<j (u=i, v=j), joint[u,v,a,b] is count(x=a, y=b)
                    # if i>j (u=j, v=i), joint[u,v,b,a] is count(x=b, y=a) -> Same thing symmetric
                    if i < j:
                        c_xy = joint[u, v, a, b]
                    else:
                        c_xy = joint[u, v, b, a]

                    if c_xy > 0:
                        # O(1) Lookup
                        mi_val += lut[c_xy, c_x, c_y]

            row_mis[j] = mi_val

        # 3. Ratio Logic (Inline Sort)
        row_mis = np.sort(row_mis)

        # Top Mean
        take_cnt = min(k_top, n_var)
        sum_top = 0.0
        if take_cnt > 0:
            for k in range(take_cnt):
                sum_top += row_mis[n_var - 1 - k]
        top_mean = sum_top / k_top

        # Bottom Mean
        sum_bot = 0.0
        if k_bottom > n_zeros_implicit:
            needed = k_bottom - n_zeros_implicit
            take_bot = min(needed, n_var)
            for k in range(take_bot):
                sum_bot += row_mis[k]
            bottom_mean = sum_bot / k_bottom
        else:
            bottom_mean = 0.0

        ratio = (top_mean + damping) / (bottom_mean + damping)
        if ratio > max_ratio:
            max_ratio = ratio

    return max_ratio


# ==========================================
# 3. PARALLEL SIMULATION DRIVER
# ==========================================
@njit(parallel=True, fastmath=True)
def run_parallel_sims(seed_int, n_sims, n_members, error_rate, valid_muts, lut, top_pct, bot_pct):
    scores = np.zeros(n_sims, dtype=np.float32)
    L = len(seed_int)

    # Parallel Loop: Each thread gets its own simulation
    for s in prange(n_sims):
        # 1. Local Allocation (Cheap in Numba)
        msa = np.empty((n_members, L), dtype=np.int8)
        var_inds = np.empty(L, dtype=np.int32)

        # 2. Mutate (Manual Loop for speed)
        # We fill msa and track variants in one go if possible,
        # but 2-pass is cleaner and safer.

        # Copy Seed
        for r in range(n_members):
            msa[r, :] = seed_int

        # Apply Mutations
        # We iterate cells.
        for r in range(n_members):
            for c in range(L):
                # Fast RNG
                if np.random.random() < error_rate:
                    # random choice 0-3
                    ri = np.random.randint(0, 4)
                    msa[r, c] = valid_muts[ri]

        # 3. Find Variants
        n_var = 0
        for c in range(L):
            v0 = msa[0, c]
            is_var = False
            for r in range(1, n_members):
                if msa[r, c] != v0:
                    is_var = True
                    break
            if is_var:
                var_inds[n_var] = c
                n_var += 1

        if n_var < 2:
            scores[s] = 0.0
        else:
            scores[s] = _fast_score_kernel(msa, n_members, L, var_inds, n_var, lut, top_pct, bot_pct)

    return scores


def fast_spoa(sequences, algorithm=1, mode='msa', max_initial_reads=50):
    """
    Corrected wrapper for the functional spoa.poa() binding.
    """
    # Ensure input is a standard list of strings
    if not isinstance(sequences, list):
        sequences = list(sequences)

    if not sequences:
        return "", []

    # --- PATH 1: Consensus Only (Fast / Approximate) ---
    # Use this for 'SPOA_alignment_initial' to break the 10s bottleneck.
    # We subsample to keep it fast.
    if mode == 'consensus' and len(sequences) > max_initial_reads:
        subset = random.sample(sequences, max_initial_reads)
        try:
            # Run SPOA on just the small subset
            consensus, _ = spoa.poa(subset, algorithm=algorithm)
            return consensus, []  # Return empty MSA because subset MSA won't match full input
        except Exception as e:
            print(f"SPOA Subset Error: {e}")
            return "", []

    # --- PATH 2: Full MSA (Required for MI Check) ---
    # Use this for 'check_cluster_integrity'. We CANNOT subsample here
    # because we need the MSA to match the input rows 1-to-1.
    try:
        # The library calculates everything in one go
        consensus, msa = spoa.poa(sequences, algorithm=algorithm)
        return consensus, msa
    except Exception as e:
        print(f"SPOA Full Error: {e}")
        return "", []


# ==========================================
# 4. MAIN FUNCTION
# ==========================================


# ==========================================
# 1. LOOKUP TABLE (Unchanged)
# ==========================================
@njit
def create_mi_lookup(n_members):
    lut = np.zeros((n_members + 1, n_members + 1, n_members + 1), dtype=np.float32)
    N = float(n_members)
    epsilon = 1e-12

    for c_xy in range(1, n_members + 1):
        for c_x in range(1, n_members + 1):
            for c_y in range(1, n_members + 1):
                # Valid counts only
                if c_xy <= c_x and c_xy <= c_y:
                    p_xy = c_xy / N
                    p_x = c_x / N
                    p_y = c_y / N

                    term = p_xy * np.log(p_xy / (p_x * p_y))
                    lut[c_xy, c_x, c_y] = term
    return lut


# ==========================================
# 2. ROBUST POPCOUNT
# ==========================================
@njit(inline='always')
def popcount64(x):
    """
    SWAR algorithm for 64-bit population count.
    Inline always to ensure it compiles to optimal assembly.
    """
    x = np.uint64(x)
    m1 = np.uint64(0x5555555555555555)
    m2 = np.uint64(0x3333333333333333)
    m4 = np.uint64(0x0f0f0f0f0f0f0f0f)

    x -= (x >> np.uint64(1)) & m1
    x = (x & m2) + ((x >> np.uint64(2)) & m2)
    x = (x + (x >> np.uint64(4))) & m4
    x = (x * np.uint64(0x0101010101010101)) >> np.uint64(56)
    return int(x)


# ==========================================
# 3. BITWISE SCORE KERNEL (Memory Safe)
# ==========================================
@njit(fastmath=True)
def _score_bitwise_internal(bitmaps, n_chunks, L, n_var, lut, top_pct, bottom_pct):
    """
    Computes MI from pre-filled bitmaps.
    bitmaps shape: (4, n_var, n_chunks) OR (4, L, n_chunks)
    We only iterate up to n_var.
    """
    # 1. Marginals (Count 1s per column per base)
    # Stack-allocated small array
    marginals = np.zeros((n_var, 4), dtype=np.int32)

    for i in range(n_var):
        for b in range(4):
            count = 0
            for k in range(n_chunks):
                val = bitmaps[b, i, k]
                if val > 0:
                    count += popcount64(val)
            marginals[i, b] = count

    # 2. Pairwise MI
    row_mis = np.zeros(n_var, dtype=np.float32)
    max_ratio = 0.0
    damping = 0.01

    k_top = max(1, int(L * top_pct))
    k_bottom = max(1, int(L * bottom_pct))
    n_zeros_implicit = L - n_var

    for i in range(n_var):
        # We calculate row i against all j
        for j in range(n_var):
            if i == j:
                row_mis[j] = 0.0
                continue

            mi_val = 0.0

            # 4x4 Base Loop
            for b1 in range(4):
                c_x = marginals[i, b1]
                if c_x == 0: continue

                for b2 in range(4):
                    c_y = marginals[j, b2]
                    if c_y == 0: continue

                    # Bitwise Intersection
                    c_xy = 0
                    for k in range(n_chunks):
                        bits = bitmaps[b1, i, k] & bitmaps[b2, j, k]
                        if bits > 0:
                            c_xy += popcount64(bits)

                    if c_xy > 0:
                        mi_val += lut[c_xy, c_x, c_y]

            row_mis[j] = mi_val

        # 3. Ratio Logic (Inline)
        # Note: sorting a small array (20-50 floats) is very fast
        row_mis = np.sort(row_mis)

        # Top Mean
        take_top = min(k_top, n_var)
        sum_top = 0.0
        if take_top > 0:
            for k in range(take_top):
                sum_top += row_mis[n_var - 1 - k]
        top_mean = sum_top / k_top

        # Bottom Mean
        sum_bot = 0.0
        if k_bottom > n_zeros_implicit:
            needed = k_bottom - n_zeros_implicit
            take_bot = min(needed, n_var)
            for k in range(take_bot):
                sum_bot += row_mis[k]
            bottom_mean = sum_bot / k_bottom
        else:
            bottom_mean = 0.0

        ratio = (top_mean + damping) / (bottom_mean + damping)
        if ratio > max_ratio:
            max_ratio = ratio

    return max_ratio


# ==========================================
# 4. PARALLEL SIMULATIONS (Crash Proof)
# ==========================================
@njit(parallel=True, fastmath=True)
def run_bitwise_sims_robust(seed_int, n_sims, n_members, error_rate, valid_muts, lut, top_pct, bot_pct):
    scores = np.zeros(n_sims, dtype=np.float32)
    L = len(seed_int)

    # Calculate chunks once
    n_chunks = (n_members + 63) // 64

    # Parallel Loop
    for s in prange(n_sims):
        # ALLOCATIONS:
        # We allocate MAX size buffers (based on L) to avoid dynamic resizing.
        # This stability prevents the allocator crash.

        # MSA Buffer: (N, L)
        msa = np.empty((n_members, L), dtype=np.int8)

        # Variant Index Buffer: (L)
        var_inds = np.empty(L, dtype=np.int32)

        # Bitmap Buffer: (4, L, n_chunks) - Fixed size!
        # Initialized to zero
        bitmaps = np.zeros((4, L, n_chunks), dtype=np.uint64)

        # 1. Reset & Mutate
        # Copy Seed
        for r in range(n_members):
            msa[r, :] = seed_int

        # Apply Random Mutations
        for r in range(n_members):
            for c in range(L):
                if np.random.random() < error_rate:
                    ri = np.random.randint(0, 4)
                    msa[r, c] = valid_muts[ri]

        # 2. Identify Variants & Pack Bits Inline
        # We combine these steps to maximize cache usage
        n_var = 0

        for c in range(L):
            # Check variance
            v0 = msa[0, c]
            is_var = False
            for r in range(1, n_members):
                if msa[r, c] != v0:
                    is_var = True
                    break

            if is_var:
                # Store this column index in our map (for logic consistency)
                var_inds[n_var] = c

                # Pack this column into bitmaps immediately
                # n_var is the index in the bitmap array
                for r in range(n_members):
                    char = msa[r, c]
                    if char > 3: continue  # Skip gaps/N

                    chunk = r // 64
                    bit = r % 64

                    # Set bit
                    bitmaps[char, n_var, chunk] |= (np.uint64(1) << np.uint64(bit))

                n_var += 1

        # 3. Score
        if n_var < 2:
            scores[s] = 0.0
        else:
            scores[s] = _score_bitwise_internal(bitmaps, n_chunks, L, n_var, lut, top_pct, bot_pct)

    return scores


# ==========================================
# 5. INTEGRATION
# ==========================================
def check_cluster_integrity_mi_hyper_optimized(cluster_df, expected_error_rate, timer, verbose=0):
    timer.start("MI_Integrity_Check_Total")

    barcodes = cluster_df['Barcode'].tolist()
    inserts = cluster_df['Insert'].tolist()
    n_members = len(barcodes)

    if n_members < 10:
        timer.stop("MI_Integrity_Check_Total")
        return False

    # 1. Alignment
    timer.start("SPOA_alignment")
    _, msa_barcodes = spoa.poa(barcodes, algorithm=1)
    _, msa_inserts = fast_spoa(inserts, algorithm=1)
    timer.stop("SPOA_alignment")

    msa_strings = [b + "-" * 1 + i for b, i in zip(msa_barcodes, msa_inserts)]
    msa_int, vocab_size = encode_msa(msa_strings)

    alphabet = "ACGTN-"
    char_to_int = {c: i for i, c in enumerate(alphabet)}
    valid_muts = np.array([char_to_int[c] for c in "ACGT"], dtype=np.int8)

    # 2. LUT
    lut = create_mi_lookup(n_members)

    # 3. Real Score
    timer.start("Real_Score")
    msa_reduced, kept_indices = remove_invariant_columns(msa_int)
    n_var_real = len(kept_indices)
    L_real = msa_int.shape[1]

    if n_var_real < 2:
        real_score = 0.0
    else:
        # Pack real data
        n_chunks = (n_members + 63) // 64
        # Allocate fixed L size or just n_var_real size (single thread is safe)
        bitmaps = np.zeros((4, n_var_real, n_chunks), dtype=np.uint64)
        msa_contig = np.ascontiguousarray(msa_reduced)

        # Manual pack for real data
        for i in range(n_var_real):
            # msa_reduced already contains only variant columns
            # so we iterate 0..n_var_real directly
            for r in range(n_members):
                char = msa_contig[r, i]
                if char > 3: continue
                chunk = r // 64
                bit = r % 64
                bitmaps[char, i, chunk] |= (np.uint64(1) << np.uint64(bit))

        real_score = _score_bitwise_internal(bitmaps, n_chunks, L_real, n_var_real, lut, 0.01, 0.50)
    timer.stop("Real_Score")

    # 4. Simulations
    timer.start("Sim_Parallel")
    seed_idx = np.random.randint(n_members)
    seed_int = np.ascontiguousarray(msa_int[seed_idx])

    sim_scores = run_bitwise_sims_robust(
        seed_int,
        100,
        n_members,
        expected_error_rate,
        valid_muts,
        lut,
        0.05,
        0.50
    )
    timer.stop("Sim_Parallel")

    threshold_99 = np.percentile(sim_scores, 99)
    is_outlier = real_score > threshold_99

    if is_outlier and verbose >= 1:
        print(f"!!! Outlier: Real={real_score:.4f} > P99={threshold_99:.4f}")

    timer.stop("MI_Integrity_Check_Total")
    return is_outlier


def calculate_msa_mi_and_top_cols(barcodes, inserts, verbose=0, mi_exclusion_distance=3):
    timer.start("SPOA_alignment_initial")
    # FIX: Use algorithm=1 (Global) to match the integrity check logic.
    _, msa_barcodes = fast_spoa(barcodes, algorithm=1, mode='msa')
    _, msa_inserts = fast_spoa(inserts, algorithm=1, mode='msa')
    timer.stop("SPOA_alignment_initial")

    msa_strings = [b + "-" * 1 + i for b, i in zip(msa_barcodes, msa_inserts)]
    len_msa_inserts = len(msa_inserts[0])
    msa_int, vocab_size = encode_msa(msa_strings)

    timer.start("Filter_Invariant")
    msa_reduced, kept_indices = remove_invariant_columns(msa_int)
    timer.stop("Filter_Invariant")

    if msa_reduced.shape[1] < 2:
        return msa_int, np.zeros((0, 0)), np.array([]), len_msa_inserts, kept_indices

    one_hot, P_i = get_marginals(msa_reduced, vocab_size)
    mi_matrix_reduced = compute_mi_scores_optimized_v2(one_hot, P_i, timer)

    # ==============================================================================
    # [NEW] MASK LOCAL NEIGHBORS (Band Masking)
    # ==============================================================================
    # We must mask interactions between columns that are physically close in the MSA
    # to prevent local sequencing errors (indels/homopolymers) from dominating MI.

    # 1. Compute physical distance matrix using broadcasting on kept_indices
    #    kept_indices is shape (N_cols,), this creates (N_cols, N_cols)
    physical_dist_matrix = np.abs(kept_indices[:, None] - kept_indices[None, :])

    # 2. Create boolean mask: True where distance is smaller than threshold
    #    This automatically covers the diagonal (dist=0)
    local_noise_mask = physical_dist_matrix < mi_exclusion_distance

    # 3. Apply mask: Set MI to 0.0 for neighbors
    mi_matrix_reduced[local_noise_mask] = 0.0
    # ==============================================================================

    # Note: We still pass mi_exclusion_distance to get_elbow_columns if it needs it
    # for other logic, but the matrix is now clean.
    top_cols_reduced = get_elbow_columns(mi_matrix_reduced, mi_exclusion_distance, verbose=verbose)

    top_cols_original = kept_indices[top_cols_reduced] if len(top_cols_reduced) > 0 else np.array([])
    return msa_int, mi_matrix_reduced, top_cols_original, len_msa_inserts, kept_indices


def cluster_msa_subset(msa_subset, expected_error_rate, error_rate_multiplier=None, jump_thresh=None):
    timer.start("hamming_distance_matrix")
    dist_matrix_subset = pairwise_distances(msa_subset, metric='hamming')
    timer.stop("hamming_distance_matrix")

    if dist_matrix_subset.shape[0] == 1:
        return dist_matrix_subset, np.array([0])

    if error_rate_multiplier is not None:
        d_thresh = expected_error_rate * error_rate_multiplier
    else:
        temp = dist_matrix_subset.copy()
        np.fill_diagonal(temp, np.inf)
        d_thresh = np.min(temp) * jump_thresh + expected_error_rate

    timer.start("agglomerative_clustering")
    agg = AgglomerativeClustering(metric="precomputed", linkage="single", distance_threshold=d_thresh, n_clusters=None)
    labels = agg.fit_predict(dist_matrix_subset)
    timer.stop("agglomerative_clustering")
    return dist_matrix_subset, labels


def cluster_barcode_insert_pairs(barcodes, inserts, mi_exclusion_distance, expected_error_rate,
                                 error_rate_multiplier=None, jump_thresh=None, verbose=0, subset_msa=True):
    msa_int, mi_matrix, top_cols, len_msa_inserts, _ = calculate_msa_mi_and_top_cols(barcodes, inserts, verbose,
                                                                                     mi_exclusion_distance)
    msa_subset = msa_int[:, top_cols] if (subset_msa and len(top_cols) > 0) else msa_int
    dist_matrix_subset, labels = cluster_msa_subset(msa_subset, expected_error_rate, error_rate_multiplier, jump_thresh)
    return msa_int, mi_matrix, msa_subset, dist_matrix_subset, labels, len_msa_inserts


def recursive_outlier_removal(all_df, cluster_df, filtered_dist_matrix, expected_error_rate, percentile_th,
                              verbose, error_rate_multiplier=1.2, n_decoys=50):
    if len(cluster_df) < 2: return cluster_df
    m = filtered_dist_matrix.astype(float);
    m[np.triu_indices_from(m)] = np.nan
    if np.nanmax(m) <= expected_error_rate * error_rate_multiplier: return cluster_df

    agg = AgglomerativeClustering(metric="precomputed", linkage="single", distance_threshold=np.nanmedian(m),
                                  n_clusters=None)
    labels = agg.fit_predict(filtered_dist_matrix)

    largest_idx = np.where(labels == np.bincount(labels).argmax())[0]
    largest_sub_dist = filtered_dist_matrix[np.ix_(largest_idx, largest_idx)]
    best_core = cluster_df['Insert'].iloc[largest_idx[np.argmin(np.mean(largest_sub_dist, axis=1))]]

    timer.start("Fuzzy_CDIST_Calculations")
    all_ratios = process.cdist([best_core], list(cluster_df['Insert']), scorer=fuzz.ratio, dtype=np.float32)[0]
    timer.stop("Fuzzy_CDIST_Calculations")

    outlier_pos = np.argmin(all_ratios)

    # Decoy Check
    timer.start("Fuzzy_CDIST_Calculations")
    decoy_ratios = \
        process.cdist([cluster_df['Insert'].iloc[outlier_pos]], all_df['Insert'].sample(n=n_decoys), scorer=fuzz.ratio,
                      dtype=np.float32)[0]
    timer.stop("Fuzzy_CDIST_Calculations")

    if all_ratios[outlier_pos] >= np.percentile(decoy_ratios, percentile_th):
        return cluster_df

    keep_mask = np.arange(len(filtered_dist_matrix)) != outlier_pos
    return recursive_outlier_removal(all_df, cluster_df.iloc[keep_mask],
                                     filtered_dist_matrix[np.ix_(keep_mask, keep_mask)],
                                     expected_error_rate, percentile_th, verbose, error_rate_multiplier, n_decoys)


# ==========================================
# RECURSIVE MI SPLITTING (NEW)
# ==========================================

def recursive_mi_splitting(cluster_df, all_df, expected_error_rate, timer, verbose, depth=0, max_depth=3):
    """
    Recursively splits a cluster if MI check indicates it is a mix of sequences.
    Splitting is based on clustering the 'top columns' (highest MI interactions).
    """
    # 1. Fast Integrity Check
    is_suspicious = check_cluster_integrity_mi_hyper_optimized(cluster_df, expected_error_rate, timer, verbose=verbose)

    # Base cases: Clean cluster, max depth, or too small to split
    if not is_suspicious:
        cluster_df = cluster_df.copy()
        cluster_df['mi_warning'] = False
        return [cluster_df]

    if depth >= max_depth:
        if verbose >= 1: print(f"  > Max depth reached ({depth}). Returning suspicious cluster.")
        cluster_df = cluster_df.copy()
        cluster_df['mi_warning'] = True
        return [cluster_df]

    if len(cluster_df) < 6:
        if verbose >= 1: print(f"  > Cluster too small to split ({len(cluster_df)}). Returning.")
        cluster_df = cluster_df.copy()
        cluster_df['mi_warning'] = True
        return [cluster_df]

    if verbose >= 1:
        print(f"  > Splitting cluster (Depth {depth}) - Size {len(cluster_df)}")

    # 2. Identify Important Columns for Splitting
    barcodes = cluster_df['Barcode'].tolist()
    inserts = cluster_df['Insert'].tolist()

    try:
        # Use existing function to get top columns
        msa_int, mi_matrix, top_cols, _, kept_indices = calculate_msa_mi_and_top_cols(barcodes, inserts, verbose=0,
                                                                                      mi_exclusion_distance=3)
    except Exception as e:
        print(f"Error in MI split calc: {e}")
        cluster_df = cluster_df.copy()
        cluster_df['mi_warning'] = True
        return [cluster_df]

    if verbose >= 2:
        print(f"    > Initial Top Cols found: {len(top_cols)}")

    # 3. Split based on Top Columns
    if len(top_cols) < 1:
        # Fallback 1: Try top 3 MI columns
        if mi_matrix is not None and mi_matrix.size > 0:
            if verbose >= 1: print("    > Elbow failed, forcing fallback split on Top 3 MI columns.")
            mi_sums = mi_matrix.sum(axis=1)
            fallback_indices = np.argsort(mi_sums)[-3:]
            top_cols = kept_indices[fallback_indices]

        # Fallback 2: If still no top cols, use ALL variant columns.
        # If MI is high, the signal MUST be in the variants somewhere.
        if len(top_cols) < 1 and len(kept_indices) > 0:
            if verbose >= 1: print("    > Fallback failed. Using ALL variant columns for split.")
            top_cols = kept_indices

        # If still nothing (no variants at all?), abort.
        if len(top_cols) < 1:
            if verbose >= 1: print("    > No variant columns found (identical sequences?). Aborting split.")
            cluster_df = cluster_df.copy()
            cluster_df['mi_warning'] = True
            return [cluster_df]

    # Extract only the high-MI columns
    msa_subset = msa_int[:, top_cols]

    # === FIX: ROBUST VECTORIZED ONE-HOT ENCODING ===
    # Use int32 to prevent any overflow issues and vectorize for speed.
    N_sub, C_sub = msa_subset.shape
    vocab_size = 6  # ACGTN-

    # Flatten for advanced indexing
    flat_vals = msa_subset.flatten()

    # Generate grid indices
    rows = np.repeat(np.arange(N_sub), C_sub)
    cols = np.tile(np.arange(C_sub), N_sub)

    # Mask invalid values to be safe
    mask = (flat_vals >= 0) & (flat_vals < vocab_size)

    valid_rows = rows[mask]
    # Calculate flattened column index for One-Hot matrix
    valid_cols = cols[mask] * vocab_size + flat_vals[mask]

    # Initialize with int32 (safe)
    one_hot = np.zeros((N_sub, C_sub * vocab_size), dtype=np.int32)
    one_hot[valid_rows, valid_cols] = 1

    timer.start("ward_splitting")
    # Ward linkage ensures we find the structural split, not just peel off an outlier
    agg = AgglomerativeClustering(n_clusters=2, linkage='ward')
    labels = agg.fit_predict(one_hot)
    timer.stop("ward_splitting")

    df_0 = cluster_df.iloc[labels == 0]
    df_1 = cluster_df.iloc[labels == 1]

    # Sanity check: Ensure neither side is strictly empty.
    # We accept 1 vs Rest splits now.
    if len(df_0) < 1 or len(df_1) < 1:
        if verbose >= 1: print(f"    > Split failed (Empty result: {len(df_0)} vs {len(df_1)}). Aborting.")
        cluster_df = cluster_df.copy()
        cluster_df['mi_warning'] = True
        return [cluster_df]

    if verbose >= 1:
        print(f"    > Split successful: {len(df_0)} vs {len(df_1)}")

    # 4. Recurse on children
    results_0 = recursive_mi_splitting(df_0, all_df, expected_error_rate, timer, verbose, depth + 1, max_depth)
    results_1 = recursive_mi_splitting(df_1, all_df, expected_error_rate, timer, verbose, depth + 1, max_depth)

    return results_0 + results_1


def full_analysis(sub_df, all_df, barcode_target, percentile_th, expected_error_rate,
                  verbose, mi_exclusion_distance, cluster_again_total_mean_thresh=2,
                  cluster_again_single_mean_thresh=3):
    if len(sub_df) < 2:
        if 'cluster' not in sub_df.columns: sub_df = sub_df.copy(); sub_df['cluster'] = -1
        return sub_df

    # Initial Clustering
    barcodes = sub_df['Barcode'].tolist()
    inserts = sub_df['Insert'].tolist()

    _, _, _, dist_matrix_subset, labels, _ = cluster_barcode_insert_pairs(
        barcodes, inserts, mi_exclusion_distance, expected_error_rate,
        error_rate_multiplier=3, verbose=0, subset_msa=False
    )

    # Unique cluster IDs to prevent collision in later steps
    labels = labels + hash(frozenset(sub_df.index)) % 10_000_000
    sub_df['cluster'] = labels

    final_clusters = []
    unique_labels = list(set(labels))

    for label in unique_labels:
        this_cluster = sub_df[sub_df['cluster'] == label]

        # 1. Outlier Removal (Existing logic)
        indices = np.where(labels == label)[0]
        filtered_dist_matrix = dist_matrix_subset[np.ix_(indices, indices)]
        processed_cluster = recursive_outlier_removal(
            all_df, this_cluster, filtered_dist_matrix,
            expected_error_rate, percentile_th, verbose
        )

        # 2. Recursive MI Check & Splitting (New Logic)
        if len(processed_cluster) >= 8 and label != -1:
            # This function returns a LIST of dataframes (1 if clean, >1 if split)
            split_results = recursive_mi_splitting(
                processed_cluster, all_df, expected_error_rate, timer, verbose
            )

            # Assign new sub-cluster IDs if split occurred
            for i, res_df in enumerate(split_results):
                if len(split_results) > 1:
                    # Append suffix to label to track lineage (e.g., 101_0, 101_1)
                    res_df = res_df.copy()
                    res_df['cluster'] = label * 100 + i  # Simple offset for sub-clusters
                final_clusters.append(res_df)

        else:
            processed_cluster = processed_cluster.copy()
            processed_cluster['mi_warning'] = False
            final_clusters.append(processed_cluster)

    if not final_clusters:
        return pd.DataFrame(columns=sub_df.columns)

    result = pd.concat(final_clusters)
    timer.report()
    return result


def filter_result_df(result, barcode_target, verbose,
                     filter_target=True,
                     min_barcode_ratio_hq=0.975):
    result = result.copy()
    unique_clusters = result['cluster'].unique()
    consensus_map = {}

    try:
        timer.start("filter_consensus_generation")
    except:
        pass

    for cl in unique_clusters:
        subset = result.loc[result['cluster'] == cl, 'Barcode']
        consensus_map[cl] = get_poa_consensus(subset)

    try:
        timer.stop("filter_consensus_generation")
    except:
        pass

    result['cluster_consensus_bc'] = result['cluster'].map(consensus_map)

    if filter_target:
        result = result[
            (result['cluster_consensus_bc'] == barcode_target) |
            ((result.groupby('cluster')['Barcode'].transform('size') == 2) &
             (result.groupby('cluster')['Barcode'].transform(lambda x: (x == barcode_target).any()))
             )
            ]

    if result.empty: return result

    def get_mean_consistency(df_chunk):
        consensus = df_chunk['cluster_consensus_bc'].iloc[0]
        if not consensus: return 0.0
        scores = [fuzz.ratio(b, consensus) for b in df_chunk['Barcode']]
        return np.mean(scores) / 100.0

    cluster_scores = result.groupby('cluster').apply(get_mean_consistency, include_groups=False)
    result['cluster_mean_barcode_sim'] = result['cluster'].map(cluster_scores)
    result['barcodes_consistent'] = result['cluster_mean_barcode_sim'] >= min_barcode_ratio_hq

    result['cluster_consensus_insert'] = result.groupby('cluster')['Insert'].transform(get_poa_consensus)

    unique = [u for u in sorted(result['cluster'].unique()) if u != -1]
    mapping = {old: new for new, old in enumerate(unique)}
    result['cluster'] = result['cluster'].map(lambda x: mapping.get(x, -1))

    return result


if __name__ == "__main__":
    # ==== PARAMS =====

    verbose = 1
    percentile_th = 95
    expected_error_rate = 0.03
    PLOTS = False
    mi_exclusion_distance = 3
    n_sims = 150
    sim_percentile = 0.1

    # ==== RUN =====

    barcode_counts = df['Barcode'].value_counts()
    all_barcodes = np.array(barcode_counts.index.tolist())

    for i, (top_bc, count) in enumerate(barcode_counts.head(20).items()):

        if verbose >= 1:
            print(f"\n\n=============\n{top_bc}")
            print(f"Processing Top Barcode #{i + 1}: {top_bc} (Count: {count})")

        # top_bc = "TATCATCC" # Force specific barcode if needed for debug

        # RapidFuzz
        scores = process.cdist([top_bc], all_barcodes, scorer=fuzz.ratio, dtype=np.float32)[0]
        candidate_bcs = all_barcodes[np.where(scores > 88)[0]]
        filtered_df = df[df['Barcode'].isin(candidate_bcs)].copy()

        print(f"length of filtered df: {len(filtered_df)}")

        # filtered_df = pd.read_csv("~/Downloads/fml.csv")

        if filtered_df.empty: continue

        if verbose >= 2:
            print(f"  > RapidFuzz gathered {len(filtered_df)} reads")

        result_df = full_analysis(filtered_df, df, top_bc,
                                  percentile_th=percentile_th, verbose=verbose, expected_error_rate=expected_error_rate,
                                  mi_exclusion_distance=mi_exclusion_distance)

        print(result_df)

        result_df2 = filter_result_df(result_df, top_bc, verbose)

        if verbose >= 1:
            # Create the count DataFrame first (using the recommended groupby method)
            combo_counts_df = result_df2.groupby(
                ['ID', 'cluster', 'mi_warning', 'barcodes_consistent']).size().reset_index(name='count')

            # Now, sort the DataFrame
            combo_counts_df.sort_values(by=['cluster', 'count'], ascending=[True, False], inplace=True)

            # Print the final result
            print(combo_counts_df)

    result_df2
