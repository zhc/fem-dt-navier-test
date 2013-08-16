from dolfin import *
import math
import sys
if len(sys.argv) < 5:
	print "Usage: navier.py nonlin nx dt reynolds"
	exit(0)

nx = 30
dt = 0.01
eps = 0.0001
T = 10
Reynolds = 500
nonlin = 1

nonlin = int(sys.argv[1])
nx = int(sys.argv[2])
dt = float(sys.argv[3])
Reynolds = int(sys.argv[4])
print "Using", nonlin, nx, dt, Reynolds

parameters["std_out_all_processes"] = False;
mesh = UnitSquareMesh(nx, nx)
U = VectorFunctionSpace(mesh, "CG", 2)
P = FunctionSpace(mesh, "CG", 1)
W = U*P
def domain_top(x, on_boundary):
    return on_boundary and x[1] > 1 - DOLFIN_EPS  
def domain_walls(x, on_boundary):
    return on_boundary and x[1] < 1 - DOLFIN_EPS    
bc0 = DirichletBC(W.sub(0), Constant((0,0)), domain_walls)
bc1 = DirichletBC(W.sub(0), Constant((1,0)), domain_top)
bcs = [bc0, bc1]    
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
uk = Function(U)
u0 = Function(U)
f = Constant((0, 0))
tau = Constant(dt)
Re = Constant(Reynolds)
a = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(grad(p), v)*dx + div(u)*q*dx + inner(grad(uk)*u, v)*dx 
L = inner(f, v)*dx + 1/tau*inner(u0, v)*dx
if nonlin==1:
	#picard 1
	pass
if nonlin==2:
	#picard 2
	a = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(grad(p), v)*dx + div(u)*q*dx + inner(grad(u)*uk, v)*dx 
	L = inner(f, v)*dx + 1/tau*inner(u0, v)*dx
elif nonlin==3:
	#picard 3
	a = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(grad(p), v)*dx + div(u)*q*dx
	L = -inner(grad(uk)*uk, v)*dx + inner(f, v)*dx + 1/tau*inner(u0, v)*dx
elif nonlin==4:
	#newton
	a = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(grad(p), v)*dx + div(u)*q*dx + inner(grad(u)*uk, v)*dx + inner(grad(uk)*u, v)*dx
	L = inner(f, v)*dx + inner(grad(uk)*uk, v)*dx + 1/tau*inner(u0, v)*dx
else:
	print "Invalid nonlin parameter:", nonlin
	exit(-1)
U = Function(W)
fvelo = File("result_%s_%d_%f_%d/velocity.pvd" % (nonlin, nx, dt, Reynolds))
#fpress = File("result_%s_%d_%f/pressure.pvd" % (nonlin, nx, dt))
t = 0
velo, press = U.split()
velo_err = 10
while t <= T:
	print "t = %lf" % t
	uk.assign(u0)
	i = 0
	velo_err = 10
	while velo_err > eps:
		solve(a == L, U, bcs=bcs)
		(velo, press) = U.split()
		M = inner((uk - velo),(uk - velo))*dx
		velo_err = assemble(M, mesh=mesh)
		print "Nonlinear convergence =", velo_err
		i = i + 1
		uk.assign(velo)
	fvelo << velo
	M = inner((u0 - velo),(u0 - velo))*dx
	velo_err = assemble(M, mesh=mesh)
	print "Stationary convergence =", velo_err
	u0.assign(velo)
	t += dt