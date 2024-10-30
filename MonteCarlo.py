# Monte Carlo method of Estimating $\ln(2)$
# Ruize Li, rl737, Jesus College, Cambridge

# The area under the function y = 1, x from 1 to 2 is $1$, and the area
# under the function y = 1/x, x from 1 to 2 is $\ln(2)$. Therefore, the
# ratio of the areas is $\ln(2)$. We can use the Monte Carlo method to
# estimate this ratio. 

# Importing necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import random

# Setting the number of iterations
N = 10000
print("Parameter: N =", str(N))

# Array of iterations and results
iterations = []
results = []

# Initializing count_in to 0
count_in = 0

# Running the Monte Carlo method
for i in range(N):
    # Generating random values of x and y, x in [1,2] and y in [0,1]
    x = 1 + random.random()
    y = random.random()

    # Checking if y is less than or equal to 1/x
    cond = y
    outcome = 1 if cond <= 1/x else 0

    # Updating the count_in and fraction_in
    count_in += outcome
    fraction_in = count_in/(i+1)

    # Appending the result and iteration to the arrays
    results.append(fraction_in)
    iterations.append(i+1)

    # Printing the progress
    print("Location: "+str(outcome)+"\t"+str(x)+"\t"+str(y)+"\t"+str(count_in)+
          "\t"+str(i)+"\t"+str(fraction_in))

# Find out how many points are required to get an answer accurate to 1 percent
for i in range(N):
    if results[-i-1] - np.log(2) > 0.01 * np.log(2):
        print("Number of iterations required to get an answer accurate to 1 percent:", N-i+1)
        # Stop the loop
        break

# Visualizing the results
fig = plt.figure()
plt.plot(iterations,results,'k-', label='Numerical $\ln(2)$')
plt.plot([0,iterations[-1]],[np.log(2),np.log(2)],'r-', label='Analytical $\ln(2)$')
plt.grid(True)
plt.legend()
plt.title('Monte Carlo Method of $\ln(2)$')
plt.xlabel('Iteration')
plt.ylabel('Result')
plt.show()