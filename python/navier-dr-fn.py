from dolfin import *
import math
import sys
if len(sys.argv) < 5:
	print "Usage: navier-dr-div.py nx dt maxt reynolds"
	exit(0)

nx = 30
dt = 0.01
eps = 0.0000001
T = 2.5
Reynolds = 500

nx = int(sys.argv[1])
dt = float(sys.argv[2])
T = float(sys.argv[3])
Reynolds = int(sys.argv[4])
print "Using", nx, dt, T, Reynolds

parameters["std_out_all_processes"] = False;
mesh = UnitSquareMesh(nx, nx)
U = VectorFunctionSpace(mesh, "CG", 2)
P = FunctionSpace(mesh, "CG", 1)
def domain_top(x, on_boundary):
    return on_boundary and x[1] > 1 - DOLFIN_EPS  
def domain_walls(x, on_boundary):
    return on_boundary and x[1] < 1 - DOLFIN_EPS    
bc0 = DirichletBC(U, Constant((0,0)), domain_walls)
bc1 = DirichletBC(U, Constant((1,0)), domain_top)
bcs = [bc0, bc1]    
u = TrialFunction(U)
v = TestFunction(U)
p = TrialFunction(P)
q = TestFunction(P)

u1 = Function(U)
u0 = Function(U)
u12 = Function(U)
ui = Function(U)

p0 = Function(P)
p1 = Function(P)

f = Constant((0, 0))
tau = Constant(dt)
Re = Constant(Reynolds)


#divergence
# F1 = 1/tau*inner(u - u0, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(div(outer(u, u0)), v)*dx - inner(f, v)*dx + inner(grad(p0), v)*dx
#F1 = 1/tau*inner(u - u0, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(div(outer(u0, u0)), v)*dx - inner(f, v)*dx + inner(grad(p0), v)*dx
F1 = 1/tau*inner(u - u0, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(grad(u)*u0, v)*dx - inner(f, v)*dx + inner(grad(p0), v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

a2 = -inner(grad(p), grad(q))*dx
L2 = 1.0/tau*div(u12)*q*dx - inner(grad(p0), grad(q))*dx

a3 = 1/tau*inner(u, v)*dx
L3 = 1/tau*inner(u12, v)*dx - inner(grad(p1-p0), v)*dx

fvelo = File("result_%d_%f_%d/velocity.pvd" % (nx, dt, Reynolds))
t = 0
velo_err = 10
while t <= T:
	print "t = %lf" % t
	i = 0
	solve(a1 == L1, u12, bcs=bcs)
	solve(a2 == L2, p1)
	solve(a3 == L3, u1, bcs=bcs)
	fvelo << u1
	M = inner((u0 - u1),(u0 - u1))*dx
	velo_err = assemble(M, mesh=mesh)
	print "Stationary convergence =", velo_err
	u0.assign(u1)
	p0.assign(p1)
	t += dt