Linear Equations Solver (Gauss & Gauss-Jordan)

This Python program provides a complete, step-by-step solution for systems of linear equations using Gauss Elimination and Gauss-Jordan Elimination.
It supports unique solutions, no solutions, and infinitely many solutions with parametric forms.
Features:
Accepts any n × (n+1) augmented matrix from the user

Implements:
1. Gauss Elimination (Forward Elimination + Back Substitution)
2. Gauss-Jordan Elimination (RREF)
3. Partial pivoting for numerical stability
4. Detailed matrix printing after every transformation

Automatically detects:
✔ Unique solution
✗ No solution (inconsistent system)
∞ Infinitely many solutions (dependent system)

🧮 What the Program Does
-> Takes number of equations and augmented matrix input
-> Allows the user to choose between:
1. Gauss Elimination
2. Gauss-Jordan Elimination

-> Shows each row operation in a formatted matrix view
-> Computes the final solution set

Displays:
Final triangular or RREF matrix

Back-substitution steps (for Gauss)

Parametric representation (for infinite solutions)
