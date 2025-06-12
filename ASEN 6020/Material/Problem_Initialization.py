###############
#READ ME FIRST
#This file is a Python file, and it uses the PyIPOPT interface with IPOPT
#This file contains the problem initialization
#The functions for the actual problem are in a second file, it will be clear below where they are linked to the problem (typically, the problem functions will be in a low-level language to benefit the run-time)
#The commands will (of course) be slightly different if you interface with IPOPT from a different language, but ultimately the commands will be almost identical
#Similarly, if you access a different optimizer it will be different as well
#Nonetheless, hopefully this shows how few commands are required to start one of these optimizers, and what the structure of one of these implementations looks like
#If you have a Python distribution with the proper modules (such as Anaconda), and compiled the IPOPT library as well as the PyIPOPT interface, this file will execute without problem
#You probably don't have all oft he above though :-) Just read this to get an idea of how the interface works
#It might be nicer to read this with a text editor that understands Python syntax, SciTE is a small & freely available example
#Everything here should correspond with the POPSICLE problem from the slides
######################

#Import the relevant Python modules
import os, sys
import numpy
import pyipopt
import Problem_Functions as PF

#Setting this global variable to 1 uses the sparse formulation, any other value uses the full Jacobian
sparse=1


def testIP():
	
	#This just affects what PyIPOPT will print
	pyipopt.set_loglevel(0)
	
	#Set the number of variables
	nvar = 3
	#Make arrays containing the lower and upper bounds for the variables
	x_L = numpy.ones((nvar), dtype=numpy.float_) * 0.0
	x_U = numpy.ones((nvar), dtype=numpy.float_) * 5.0
	
	#Set number of constraints
	ncon = 2
	
	#Make arrays containing the lower and upper bounds for the constraints
	g_L = numpy.array([0.0, 30.0])
	g_U = numpy.array([2.5, 100.0]) 
	
	#Depending on whether sparsity is used, set the number of values in the Jacobian
	if sparse==1:
		nnzj=4
	else:
		nnzj = 6
	
	#Ignore these inputs, they're for more advanced use of IPOPT, but are expected by the interface
	pi0 = numpy.array([0.0, 0.0])
	nnzh = 10 

	#Instantiate the optimization problem using the values defined so far
	#The last four inputs are the cost function, cost derivatives, constraints, and constraint derivatives, respectively
	#These functions are in a separate file, which was imported at the top of this file
	nlp = pyipopt.create(nvar, x_L, x_U, ncon, g_L, g_U, nnzj, nnzh, PF.eval_f, PF.eval_grad_f, PF.eval_g, PF.eval_jac_g)

	
	#These are some example changes in IPOPT settings, but I'm actually just resetting them to their original value (its just a demonstration of where/how to change options)
	nlp.num_option("tol", 10**(-6)) 
	nlp.int_option("file_print_level", 5) 
	
	
	#Create some sort of initial guess
	x0 = numpy.array([0.0, 2., 2.0])
	
	
	#Run the optimizer for the defined problem, using the aboe initial guess
	print "Going to call solve"
	x, zl, zu, constraint_multipliers, obj, status = nlp.solve(x0)
	
	#Close the problem
	nlp.close()
	
	#Print a bunch of interesting values
	print "Objective value"
	print "f(x*) =", obj #This value should be (roughly): 60.7519130075
	print x0
	print x #These values should be (roughly): [  1.78191104e-05   3.79143762e+00   2.50000002e+00]
	
	
#This bottom part means that the above function will execute if you run this from the command line as follows: python Problem_Initialization.py
#This is just a Python specific thing, you can totally ignore this
if __name__ == "__main__":
    testIP()
     