# Session 1
# Ruize Li, rl737
# Jesus College Cambridge

# SimpleAdder_plus_PolynomialSolver.py consists of two parts:
# 1. A simple adder program in Python;
# 2. A program to find the zeros of a polynomial function in Python.

# Choose one task:
c = input("Which one do you want to do? (1/2) \n 1. Simple Adder \n 2. Find Zeros of Polynomial Function")
if c == "1":
    # continue with the simple adder program
    # The first part is the code for a simple adder program in Python.
    # There is no library or module needed for this program.
    # Here is the function for a simple adder:
    def add(a, b):
        return a + b

    # Input the first number to add:
    a = input("This is a simple adder that can add two NUMBERS. \n Enter the first number: ")

    # Try to convert the input to a float:
    try: 
        a = float(a)
    except ValueError:
        # Notifying the user and exiting the program:
        raise ValueError("Input is not a number")

    # Input the second number to add:
    b = input("Enter the second number: ")

    # Checking if the input is a number:
    try:
        b = float(b)
    except ValueError:
        # Notifying the user and exiting the program:
        raise ValueError("Input is not a number")

    # adding the two numbers:
    total = add(float(a), float(b))

    # printing the result:
    print("The sum of", a, "and", b, "is", total)
elif c == "2":
    # continue with the program to find the zeros of a polynomial function
    # The second part of the program is to find the zeros of a polynomial function in Python.
    # Choose to continue or exit the program:
    confirm = input("This program finds the zeros of a polynomial function. \n Do you want to continue? (y/n) ")
    if confirm.lower() != "y":
        print("Exiting the program...")
        raise SystemExit("You have chosen to exit the program.")
    else:
        # continue with the program
        pass

    # Importing the necessary libraries:
    import numpy as np

    # Determining the degree of the polynomial function:
    degree = input("Enter the degree (order) of the polynomial function: ")

    # Checking if the input is a number:
    try:
        degree = int(degree)
    except ValueError:
        # Notifying the user and exiting the program:
        raise ValueError("Input is not a number")

    # Inputting the coefficients of the polynomial function:
    coefficients = []
    for i in range(degree+1):
        coefficient = input("Enter the coefficient of x^{}: ".format(degree-i))
        # Checking if the input is a number:
        try:
            coefficient = float(coefficient)
        except ValueError:
        # Notifying the user and exiting the program:
            raise ValueError("Input is not a number")
        # Adding the coefficient to the list of coefficients:
        coefficients.append(coefficient)

    # Defining a function to find the zeros of a polynomial function:
    def find_polynomial_zeros(coefficients):
        roots = np.roots(coefficients)
        return roots

    # Finding the zeros of the polynomial function:
    zeros = find_polynomial_zeros(coefficients)

    # Checking if the zeros are real:
    roots = []
    for i in range(0,len(zeros)):
        # Choosing 1e-5 as the threshold for imaginary part:
        if abs(zeros[i].imag) < 1e-5:
            roots.append(str(round(zeros[i].real,4)))
        else:
            roots.append(str(round(zeros[i].real,4))+" + "+str(round(zeros[i].imag,4))+"i")

    # Printing the polynomial function:
    print("The polynomial function is: ")
    print(str(coefficients[0]) + " x^" + str(int(degree)) + " + " 
        + " + ".join([str(coefficients[i]) + " x^" + str(int(degree-i)) for i in range(1,degree+1)]))

    # Print the zeros of the polynomial:
    print("The zeros of the polynomial are:")
    print(", ".join(roots[i] for i in range(len(roots))))
else:
    print("Invalid input. Exiting the program...")
    raise SystemExit("You have exited the program.")

# End of the program.
print("Thank you for using the program. See you next time!")