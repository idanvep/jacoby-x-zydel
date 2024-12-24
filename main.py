def is_dominant_diagonally(A):
    r = len(A)
    for i in range(r):
        diag_elem = abs(A[i][i])  # Absolute value of diagonal element
        non_diag_sum = sum(abs(A[i][j]) for j in range(r) if j != i)  # Sum of non-diagonal elements in row
        if diag_elem < non_diag_sum:  # Check if diagonal element is smaller than sum of non-diagonal elements
            return False  # Return False if any row is not diagonally dominant
    return True  # Return True if all rows are diagonally dominant


# Function to decompose a matrix into its lower, upper, and diagonal components
def decompose_matrix(A):
    n = len(A)
    L = [[0] * n for _ in range(n)]
    U = [[0] * n for _ in range(n)]
    D = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if i > j:
                L[i][j] = A[i][j]
            elif i < j:
                U[i][j] = A[i][j]
            else:
                D[i][j] = A[i][j]

    return L, U, D


# Function to perform the Gauss-Seidel method using the decomposed matrices
def GS_method_with_decomposition(L, U, D, Y, X, epsilon=0.00001, max_iterations=1000):
    n = len(D)
    for iteration in range(max_iterations):
        X_new = X.copy()  # Create a copy of the current guess
        for j in range(n):
            # Calculate the summation of L and U parts
            summ_val = Y[j] - sum(L[j][k] * X_new[k] for k in range(j)) - sum(U[j][k] * X[k] for k in range(j + 1, n))
            X_new[j] = summ_val / D[j][j]

        # Calculate the maximum difference between the current and previous estimates
        max_diff = max(abs(X_new[i] - X[i]) for i in range(len(X)))
        X = X_new  # Update the guess with the new values

        print(f"Iteration {iteration + 1}: {X}")  # Print the current iteration and guess

        if max_diff < epsilon:
            print(f"Converged after {iteration + 1} iterations.")
            return X  # Return the guess if the solution has converged

    print("Did not converge within the maximum number of iterations")
    return X  # Return the final guess if it did not converge within the maximum iterations


# Function to perform the Jacobi method using the decomposed matrices
def jacobi_method(L, U, D, Y, X, epsilon=0.00001, max_iterations=1000):
    n = len(D)
    for iteration in range(max_iterations):
        X_new = X.copy()  # Create a copy of the current guess
        for j in range(n):
            # Calculate the summation of L and U parts
            summ_val = Y[j] - sum(L[j][k] * X[k] for k in range(n)) - sum(U[j][k] * X[k] for k in range(n))
            X_new[j] = summ_val / D[j][j]

        # Calculate the maximum difference between the current and previous estimates
        max_diff = max(abs(X_new[i] - X[i]) for i in range(len(X)))
        X = X_new  # Update the guess with the new values

        print(f"Iteration {iteration + 1}: {X}")  # Print the current iteration and guess

        if max_diff < epsilon:
            print(f"Converged after {iteration + 1} iterations.")
            return X  # Return the guess if the solution has converged

    print("Did not converge within the maximum number of iterations")
    return X  # Return the final guess if it did not converge within the maximum iterations


# Function to check if the matrix is diagonally dominant or can be made so by row swapping
def is_diagonally_dominant(matrixA):
    if is_dominant_diagonally(matrixA):
        return matrixA
    else:
        # Try to make it diagonally dominant by row swapping
        for i in range(len(matrixA)):
            for j in range(i + 1, len(matrixA)):
                matrixA[i], matrixA[j] = matrixA[j], matrixA[i]  # Swap rows
                if is_dominant_diagonally(matrixA):
                    return matrixA  # Return matrix if it becomes diagonally dominant
                # Swap back if not dominant
                matrixA[i], matrixA[j] = matrixA[j], matrixA[i]
    return None  # Return None if matrix cannot be made diagonally dominant


# Main execution
matrixA = [[4, 2, 0], [2, 10, 4], [0, 4, 5]]  # Coefficient matrix
vectorB = [2, 6, 5]  # Constants on the right-hand side of the equations
initial_guess = [0, 0, 0]  # Initial guess for the solution

# Decompose matrix A into L, U, and D
L, U, D = decompose_matrix(matrixA)

# Check for diagonal dominance
dominant_matrix = is_diagonally_dominant(matrixA)

if dominant_matrix:
    print("Matrix is diagonally dominant or made diagonally dominant by row swapping.")
else:
    print("Matrix is not diagonally dominant even after row swapping. Running limited iterations.")

# Get user input for the method choice
method_choice = input("Choose a method: Gauss-Seidel (1) or Jacobi (2): ").strip().lower()

if method_choice == '1':
    print("Using Gauss-Seidel Method:")
    GS_method_with_decomposition(L, U, D, vectorB, initial_guess)
elif method_choice == '2':
    print("Using Jacobi Method:")
    jacobi_method(L, U, D, vectorB, initial_guess)
else:
    print("Invalid choice. Please choose either 'Gauss-Seidel' or 'Jacobi'.")

