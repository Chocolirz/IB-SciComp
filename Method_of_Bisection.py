# Session 4: Bisection method
# Ruize Li, rl737, Jesus College Cambridge

# This is a Python implementation of solving quadratic equations using 
# the Lehmer-Schur method which uses the Bisection method for finding the square root
# instead of using the built-in math.sqrt() function.

def sqrt_bisection(N, tolerance=1e-6, max_iterations=1000):
    # This function finds the square root of a number N using the Bisection method.

    # Set the initial condition
    complexnum = False

    # Edge case for zero
    if N == 0:
        return 0.0

    # Check if the number is non-negative
    elif N < 0:
        N = -N
        complexnum = True

    # Initial interval [a, b]
    a = 0.0
    b = max(1.0, N)

    # Bisection method loop
    iterations = 0
    while iterations < max_iterations:
        m = (a + b) / 2  # Midpoint
        f_m = m * m - N  # Function value at midpoint
        
        # If we are within the desired tolerance, return the midpoint
        if abs(f_m) < tolerance and not complexnum:
            return m # Return the number directly
        elif abs(f_m) < tolerance and complexnum:
            return complex(0, round(m,4))  # Return a complex number with zero real part
        
        # Narrow the interval based on the sign of f(m)
        if f_m > 0:
            b = m  # The square root is in the lower half
        else:
            a = m  # The square root is in the upper half
        
        iterations += 1
    
    # Return the best estimate after the loop
    if complexnum:
        return complex(0, round((a + b) / 2,4))  # Return a complex number with zero real part
    else:
        return complex(round((a + b) / 2,4), 0)  # Return a complex number with zero imaginary part

def lehmer_schur(a, b, c):
    # This function solves the quadratic equation ax^2 + bx + c = 0 using the Lehmer-Schur method.
    if a == 0:
        raise ValueError("Coefficient 'a' cannot be zero for a quadratic equation.")

    # Step 1: Compute the necessary discriminants and constants
    p = b / a
    q = c / a

    # Step 2: Calculate the first step in the Lehmer-Schur iteration
    r = sqrt_bisection(p * p / 4.0 - q)  # Calculate sqrt((b^2 - 4ac) / 4a^2)

    # Step 3: Compute the two possible roots using Lehmer-Schur's algorithm
    root1 = -p / 2 + r
    root2 = -p / 2 - r

    # Step 4: Return the roots
    return (root1, root2)

# Example
def ask():
    # Ask the user for the coefficients and constant term of the quadratic equation
    if input("Do you want to solve the quadratic equation (y/n)") == "y":
        a = input("Enter the coefficient 'a': ")
        b = input("Enter the coefficient 'b': ")
        c = input("Enter the constant term 'c': ")
    else:
        quit()
    return [float(a), float(b), float(c)]

# Test the function
a, b, c = ask()
roots = lehmer_schur(a, b, c)
print(f"The roots of the equation are: {roots}")