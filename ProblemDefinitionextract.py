# This is the changed code from ProblemDefinition.py
#The inputs are
#func_derivative_calc which is the value from the supplied Jacobian
#func_derivative_num which is the value from the numerical derivative of the function
#both values are take from a point in the middle of the range
#In this extract I take relative difference and check if it is close enough to zero.



# now compare the actual derivatives with the calculated ones   
#Relative difference
difference = (func_derivative_calc - func_derivative_num)/func_derivative_calc
#deal with nan results, +- inf should be rejected
difference[numpy.isnan(difference)]= 0 
#same check
if not (numpy.round(difference, places) == 0).all():
    raise ValueError("analytical derivative matrix does not match numerical derivative matrix.\n Analytical is:\n" + str(func_derivative_calc)
                     + "\n Numerical is:\n" + str(func_derivative_num))
