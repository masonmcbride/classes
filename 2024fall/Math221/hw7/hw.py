import numpy as np
import matplotlib.pyplot as plt

def generate_sketching_matrix_F(k, m):
    Q, R = np.linalg.qr(np.random.randn(m,m))
    F = Q[:k,:] * np.sqrt(m / k)  # F is k x m
    return F

def generate_random_matrix_A(m, n):
    t = 10  # Number of singular values in each group
    num_groups = int(np.ceil(n / t))
    r = (1e-10) ** (1 / (num_groups - 1)) if num_groups > 1 else 1e-10

    # Generate singular values
    singular_values = []
    for j in range(num_groups):
        value = r ** j
        count = t if (len(singular_values) + t) <= n else n - len(singular_values)
        singular_values.extend([value] * count)
    true_sigma_i = np.array(singular_values)

    # Ensure the length of singular_values is n
    assert len(true_sigma_i) == n, "The length of singular values should be n."

    # Generate random orthogonal matrices U and V_T
    U, _ = np.linalg.qr(np.random.randn(m,n))
    V, _ = np.linalg.qr(np.random.randn(m,n))

    # Create the diagonal matrix Sigma
    Sigma = np.diag(true_sigma_i)
    
    # Construct the matrix A
    A = U @ Sigma @ V.T

    # Return A and the true singular values
    return A, U, true_sigma_i, V.T

import numpy as np
import matplotlib.pyplot as plt

def generate_sketching_matrix_F(k, m):
    # Generate a random m x m orthogonal matrix
    Q, _ = np.linalg.qr(np.random.randn(m, m))
    # Take the first k rows and scale
    F = Q[:k, :] * np.sqrt(m / k)  # F is k x m
    return F

def generate_random_matrix_A(m, n):
    t = 10  # Number of singular values in each group
    num_groups = int(np.ceil(n / t))
    r = (1e-10) ** (1 / (num_groups - 1)) if num_groups > 1 else 1e-10

    # Generate singular values
    singular_values = []
    for j in range(num_groups):
        value = r ** j
        count = t if (len(singular_values) + t) <= n else n - len(singular_values)
        singular_values.extend([value] * count)
    true_sigma_i = np.array(singular_values)

    # Ensure the length of singular_values is n
    assert len(true_sigma_i) == n, "The length of singular values should be n."

    # Generate random orthogonal matrices U and V_T
    U, _ = np.linalg.qr(np.random.randn(m, n))
    V, _ = np.linalg.qr(np.random.randn(n, n))

    # Create the diagonal matrix Sigma
    Sigma = np.diag(true_sigma_i)
    
    # Construct the matrix A
    A = U @ Sigma @ V.T

    # Return A and the true singular values
    return A, true_sigma_i

def perform_tests():
    for n in [50, 500]:
        for m in [2 * n, 10 * n]:
            # Generate A and true_sigma_i
            A, true_sigma_i = generate_random_matrix_A(m, n)
            for k in [int(0.1 * n), int(0.5 * n)]:
                # Initialize a dictionary to collect ratios per unique singular value
                unique_true_sigma_i = np.unique(true_sigma_i)
                ratio_dict = {sigma: [] for sigma in unique_true_sigma_i}

                for _ in range(20):
                    F = generate_sketching_matrix_F(k, m)
                    # Compute singular values of FA
                    sketched_sigma_i = np.linalg.svd(F @ A, compute_uv=False)
                    # Match lengths
                    min_len = min(len(sketched_sigma_i), len(true_sigma_i))
                    # Compute ratios
                    ratios_current = sketched_sigma_i[:min_len] / true_sigma_i[:min_len]
                    # Store ratios in the dictionary
                    for idx in range(min_len):
                        sigma = true_sigma_i[idx]
                        ratio = ratios_current[idx]
                        ratio_dict[sigma].append(ratio)

                # Prepare data for plotting
                data_to_plot = []
                positions = []
                for sigma in sorted(unique_true_sigma_i, reverse=True):
                    data_to_plot.append(ratio_dict[sigma])
                    positions.append(sigma)

                # Plotting
                print(f"{m=} {n=} {k=}")
                plt.figure(figsize=(10, 6))
                plt.boxplot(data_to_plot, positions=positions, vert=False)
                plt.xscale('log')
                plt.xlabel('Ratio of Sketched to True Singular Values')
                plt.yscale('log')
                plt.ylabel('True Singular Values')
                plt.title(f'Whisker Plot (n={n}, m={m}, k={k})')
                plt.tight_layout()

                plt.tight_layout()
                plt.show()

# Run the tests
perform_tests()