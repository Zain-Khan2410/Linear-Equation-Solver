def print_matrix(matrix, step_description=""):
    if step_description:
        print(f"\n{step_description}")
    print("-" * 50)
    for row in matrix:
        print("[", end=" ")
        for val in row:
            print(f"{val:8.2f}", end=" ")
        print("]")
    print("-" * 50)


def get_augmented_matrix():
    n = int(input("Enter the number of equations: "))
    
    print(f"\nEnter the augmented matrix ({n}x{n+1}):")
    print("(Enter each row with space-separated values, including the constant)")
    
    matrix = []
    for i in range(n):
        while True:
            try:
                row_input = input(f"Row {i+1}: ").strip().split()
                if len(row_input) != n + 1:
                    print(f"Error: Expected {n+1} values, got {len(row_input)}")
                    continue
                row = [float(x) for x in row_input]
                matrix.append(row)
                break
            except ValueError:
                print("Error: Please enter valid numbers")
    
    return matrix


def swap_rows(matrix, i, j):
    matrix[i], matrix[j] = matrix[j], matrix[i]


def find_pivot(matrix, col, start_row):
    n = len(matrix)
    max_val = abs(matrix[start_row][col])
    max_row = start_row
    
    for i in range(start_row + 1, n):
        if abs(matrix[i][col]) > max_val:
            max_val = abs(matrix[i][col])
            max_row = i
    
    return max_row if max_val > 1e-10 else -1


def gauss_elimination(matrix):
    n = len(matrix)
    m = len(matrix[0])
    
    print("\n" + "="*50)
    print("GAUSS ELIMINATION METHOD")
    print("="*50)
    print_matrix(matrix, "Initial Augmented Matrix:")
    
    step = 1
    
    for col in range(n):
        pivot_row = find_pivot(matrix, col, col)
        
        if pivot_row == -1:
            continue
        
        if pivot_row != col:
            swap_rows(matrix, col, pivot_row)
            print_matrix(matrix, f"Step {step}: Swap R{col+1} ↔ R{pivot_row+1}")
            step += 1
        
        pivot = matrix[col][col]
        
        if abs(pivot) < 1e-10:
            continue
        
        for i in range(col + 1, n):
            if abs(matrix[i][col]) < 1e-10:
                continue
            
            factor = matrix[i][col] / pivot
            for j in range(m):
                matrix[i][j] -= factor * matrix[col][j]
            
            print_matrix(matrix, f"Step {step}: R{i+1} = R{i+1} - ({factor:.2f}) × R{col+1}")
            step += 1
    
    return matrix


def gauss_jordan(matrix):
    n = len(matrix)
    m = len(matrix[0])
    
    print("\n" + "="*50)
    print("GAUSS-JORDAN METHOD")
    print("="*50)
    print_matrix(matrix, "Initial Augmented Matrix:")
    
    step = 1
    pivot_cols = []
    
    row = 0
    for col in range(n):
        if row >= n:
            break
        
        pivot_row = find_pivot(matrix, col, row)
        
        if pivot_row == -1:
            continue
        
        pivot_cols.append(col)
        
        if pivot_row != row:
            swap_rows(matrix, row, pivot_row)
            print_matrix(matrix, f"Step {step}: Swap R{row+1} ↔ R{pivot_row+1}")
            step += 1
        
        pivot = matrix[row][col]
        
        if abs(pivot - 1.0) > 1e-10:
            for j in range(m):
                matrix[row][j] /= pivot
            print_matrix(matrix, f"Step {step}: R{row+1} = R{row+1} / {pivot:.2f}")
            step += 1
        
        for i in range(n):
            if i == row or abs(matrix[i][col]) < 1e-10:
                continue
            
            factor = matrix[i][col]
            for j in range(m):
                matrix[i][j] -= factor * matrix[row][j]
            
            print_matrix(matrix, f"Step {step}: R{i+1} = R{i+1} - ({factor:.2f}) × R{row+1}")
            step += 1
        
        row += 1
    
    return matrix, pivot_cols


def back_substitution(matrix):
    n = len(matrix)
    m = len(matrix[0])
    
    rank = 0
    for i in range(n):
        has_nonzero = False
        for j in range(n):
            if abs(matrix[i][j]) > 1e-10:
                has_nonzero = True
                rank += 1
                break
        
        if not has_nonzero and abs(matrix[i][m-1]) > 1e-10:
            return None, "no_solution", None
    
    if rank < n:
        pivot_cols = []
        for i in range(n):
            for j in range(n):
                if abs(matrix[i][j]) > 1e-10:
                    pivot_cols.append(j)
                    break
        return None, "infinite_solutions", pivot_cols
    
    solution = [0.0] * n
    
    print("\n" + "="*50)
    print("BACK SUBSTITUTION")
    print("="*50)
    
    for i in range(n - 1, -1, -1):
        solution[i] = matrix[i][m-1]
        for j in range(i + 1, n):
            solution[i] -= matrix[i][j] * solution[j]
        solution[i] /= matrix[i][i]
        print(f"x{i+1} = {solution[i]:.2f}")
    
    return solution, "unique", None


def analyze_solution_gauss_jordan(matrix, pivot_cols):
    n = len(matrix)
    m = len(matrix[0])
    
    for i in range(n):
        all_zero = True
        for j in range(n):
            if abs(matrix[i][j]) > 1e-10:
                all_zero = False
                break
        
        if all_zero and abs(matrix[i][m-1]) > 1e-10:
            return None, "no_solution"
    
    rank = len(pivot_cols)
    
    if rank == n:
        solution = [matrix[i][m-1] for i in range(n)]
        return solution, "unique"
    else:
        return None, "infinite_solutions"


def print_parametric_form(matrix, pivot_cols, method):
    n = len(matrix)
    m = len(matrix[0])
    
    free_vars = [i for i in range(n) if i not in pivot_cols]
    
    if not free_vars:
        print("Error: No free variables identified")
        return
    
    free_var_names = [f"x{i+1}" for i in free_vars]
    param_names = [f"t{i+1}" for i in range(len(free_vars))]
    
    print(f"\nFree variables: {', '.join(free_var_names)}")
    print(f"Let {', '.join([f'{free_var_names[i]} = {param_names[i]}' for i in range(len(free_vars))])}")
    
    print("\n" + "="*50)
    print("PARAMETRIC SOLUTION:")
    print("="*50)
    
    for var_idx in range(n):
        if var_idx in free_vars:
            param_idx = free_vars.index(var_idx)
            print(f"x{var_idx+1} = {param_names[param_idx]}")
        else:
            row_idx = pivot_cols.index(var_idx)
            
            constant = matrix[row_idx][m-1]
            terms = []
            
            if abs(constant) > 1e-10:
                terms.append(f"{constant:.2f}")
            elif not free_vars:
                terms.append("0")
            
            for free_idx, free_var in enumerate(free_vars):
                coef = -matrix[row_idx][free_var]
                
                if abs(coef) > 1e-10:
                    if coef > 0:
                        if terms:
                            terms.append(f"+ {coef:.2f}{param_names[free_idx]}")
                        else:
                            terms.append(f"{coef:.2f}{param_names[free_idx]}")
                    else:
                        terms.append(f"- {abs(coef):.2f}{param_names[free_idx]}")
            
            if not terms:
                terms.append("0")
            
            print(f"x{var_idx+1} = {' '.join(terms)}")
    
    print("="*50)
    print(f"\nwhere {', '.join(param_names)} ∈ ℝ")
    print("="*50)


def print_solution(solution, solution_type, method, matrix=None, pivot_cols=None):
    print("\n" + "="*50)
    print("SOLUTION")
    print("="*50)
    
    if solution_type == "unique":
        print("\n✓ UNIQUE SOLUTION FOUND")
        print("-" * 50)
        for i, val in enumerate(solution):
            print(f"x{i+1} = {val:.2f}")
    
    elif solution_type == "no_solution":
        print("\n✗ NO SOLUTION (Inconsistent System)")
        print("-" * 50)
        print("The system has no solution because it contains")
        print("contradictory equations (e.g., 0 = non-zero constant)")
    
    elif solution_type == "infinite_solutions":
        print("\n∞ INFINITELY MANY SOLUTIONS (Dependent System)")
        print("-" * 50)
        
        if matrix is not None and pivot_cols is not None:
            print_parametric_form(matrix, pivot_cols, method)
        else:
            print("The system has infinitely many solutions.")
    
    print("="*50)


def main():
    while True:
        print("\n" + "="*50)
        print("LINEAR EQUATIONS SOLVER")
        print("="*50)
        print("Muhammad Zain Khan (24K-0838)")
        print("Muhammad Zamin (24K-0982)")
        print("Atta Hussain (24K-0797)")
        print("="*50)
        
        print("\nChoose the solution method:")
        print("1. Gauss Elimination")
        print("2. Gauss-Jordan Elimination")
        print("3. Exit")
        
        while True:
            choice = input("\nEnter your choice (1, 2, or 3): ").strip()
            if choice in ['1', '2', '3']:
                break
            print("Invalid choice. Please enter 1, 2, or 3.")
        
        if choice == '3':
            print("\nThank you for using Linear Equations Solver!")
            print("Goodbye! 👋")
            break
        
        matrix = get_augmented_matrix()
        matrix_copy = [row[:] for row in matrix]
        
        if choice == '1':
            result_matrix = gauss_elimination(matrix_copy)
            print_matrix(result_matrix, "Final Upper Triangular Matrix:")
            
            solution, solution_type, pivot_cols = back_substitution(result_matrix)
            print_solution(solution, solution_type, "gauss", result_matrix, pivot_cols)
        
        else:
            result_matrix, pivot_cols = gauss_jordan(matrix_copy)
            print_matrix(result_matrix, "Final Reduced Row Echelon Form (RREF):")
            
            solution, solution_type = analyze_solution_gauss_jordan(result_matrix, pivot_cols)
            print_solution(solution, solution_type, "gauss_jordan", result_matrix, pivot_cols)
        
        print("\n" + "="*50)
        continue_choice = input("\nSolve another system? (y/n): ").strip().lower()
        if continue_choice != 'y':
            print("\nThank you for using Linear Equations Solver!")
            print("Goodbye! 👋")
            break


if __name__ == "__main__":
    main()
