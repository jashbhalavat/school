###############################################
#If you opened this file first, please read the Problem_Initialization file first
###############################################

#Import the relevant modules
import pyipopt
from numpy import *
import __main__


#Compute the cost function, pretty straightforward
def eval_f(x, user_data = None):
	f= x[0]**2 + (x[1])**3+x[2]**2
	return f

#Compute the cost derivatives, also straightforward
def eval_grad_f(x, user_data = None):
	g_obj=zeros((3))
	g_obj[0]=2*x[0]
	g_obj[1]=3*(x[1])**2 
	g_obj[2]=2*x[2] 
	return g_obj
	
#Compute the constraints, yet again...straightforward!
def eval_g(x, user_data= None):

	g=zeros((2))
	g[0]=x[0]**3  + x[2]    
	g[1]=(x[1])**2+x[2]**3 
	
	return g

#Compute the constraint derivatives, a little more involved
def eval_jac_g(x, flag, user_data = None):
	
	#Check if the sparse formulation is desired
	#(normally you wouldn't be checking for this, and most likely just have the sparse formulation hardcoded)
	#(I just figured it might be nice to see both)
	if __main__.sparse==1:
		
		#This is IPOPT specific
		#It will call your Jacobian function once at the start of the problem with Flag=true, in which case it's asking for the sparsity structure
		#Hence, it is defined below
		#All other iterations will have Flag=false, asking for the value of the derivatives rather than their indexing
		if flag:
			#These can be in any order, as long as they are consistent throughout the formulation
			#See also IPOPT's documentation
			#For SNOPT/fmincon/... the sparsity formulation will be a little bit different, but the idea is the same
			iRow=array([0,1,1, 0])
			jCol=array([0, 1,2, 2])
			return iRow, jCol
			
		else:
			#Compute the derivatives and put them in the right spot
			jac_g=zeros((4))
			jac_g[0]=3*x[0]**2
			jac_g[1]=2*x[1]
			jac_g[2]=3*x[2]**2
			jac_g[3]= 1
			
			
			
			
			
	#This is for the dense Jacobian		
	else:
		if flag:
			#Again returning the indexing
			iRow=array([0, 1, 0, 1, 0, 1])
			jCol=array([0, 0, 1, 1, 2, 2])
			return iRow, jCol
			
		else:
			#I don't actually need to create a 2 dimensional array, but why not
			g_con1=zeros([2,3])
			g_con1[0,0]=3*x[0]**2
			g_con1[1,1]=2*x[1]
			g_con1[1,2]=3*x[2]**2
			g_con1[0,2]=1
			
			#Inevitably, the values end up in a 1-dimensional array to be returned to IPOPT
			g_con2=zeros((6))
			g_con2[0]=g_con1[0,0]
			g_con2[1]=g_con1[1,0]
			g_con2[2]=g_con1[0,1]
			g_con2[3]=g_con1[1,1]
			g_con2[4]=g_con1[0,2]
			g_con2[5]=g_con1[1,2]
			jac_g=g_con2
			
	return jac_g
		
