### This is a program to estimate the value of pi using collisions. 

'''
Theory:

If we have two masses, m and M = 10^n m, and m is placed between M and a wall. 
Assuming no friction and perfect elastic collisions, M is given an initial velocity v_M_0 = - 1 to the wall. 
(We have set the positive direction to be the direction away from the wall.)
m is initially at rest. 
The total number of collisions should equal to 10^m times pi. 

For each collision, the velocity should be renewed to:
u_m = ((m-M)/(m+M))*v_m + ((2M)/(m+M))*v_M;
u_M = ((M-m)/(m+M))*v_M + ((2m)/(m+M))*v_m.

There's no need to consider the time interval, the wall simply gives
v_m = - u_m;
v_M = u_M, 
and the process continues until it reaches the following condition:
0 < v_m < v_M.
'''

def estimate_pi(n):
    i = 0
    v_M = -1
    v_m = 0
    m = 1
    M = 100**n
    # run the loop
    while True:
        # collision between two blocks
        u_m = ((m-M)/(m+M))*v_m + ((2*M)/(m+M))*v_M
        u_M = ((M-m)/(m+M))*v_M + ((2*m)/(m+M))*v_m
        i = i + 1

        # sort out if a final collision with the wall occures
        if abs(u_m) < abs(u_M):
            if u_m < 0:
                i = i + 1
            break

        # collision between m and the wall
        v_m = - u_m
        v_M = u_M
        i = i + 1

        # apply the final condition
        if abs(v_m) < abs(v_M):
            break
    return i

# print the result
for i in range(1, 7):
    print('The total number of collisions subject to n = ', i, 'is:', estimate_pi(i))
