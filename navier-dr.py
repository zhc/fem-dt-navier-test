from dolfin import *
import math
import sys
if len(sys.argv) < 6:
	if (MPI.process_number() == 0):
		print "Usage: navier-dr.py N DT T CONV RE"
	exit(-1)

N = 30
DT = 0.01
T = 2.5
CONV = 'd1'
REYNOLDS = 500

N = int(sys.argv[1])
DT = float(sys.argv[2])
T = float(sys.argv[3])
CONV = sys.argv[4]
REYNOLDS = int(sys.argv[5])
if (MPI.process_number() == 0):
	print "Using", N, DT, T, CONV, REYNOLDS

eps = 0.0000001

mesh = UnitSquareMesh(N, N)
U = VectorFunctionSpace(mesh, "CG", 2)
P = FunctionSpace(mesh, "CG", 1)
def domain_top(x, on_boundary):
    return on_boundary and x[1] > 1 - DOLFIN_EPS  
def domain_walls(x, on_boundary):
    return on_boundary and x[1] < 1 - DOLFIN_EPS    
def domain_all_walls(x, on_boundary):
    return on_boundary
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

p0 = Function(P)
p1 = Function(P)

f = Constant((0, 0))
tau = Constant(DT)
Re = Constant(REYNOLDS)

a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(div(outer(u0, u)), v)*dx 
L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx
if CONV == 'd1':
	a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(div(outer(u0, u)), v)*dx 
	L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx
elif CONV == 'd2':
	a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(div(outer(u, u0)), v)*dx 
	L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx
elif CONV == 'd3':
	a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx 
	L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx - inner(div(outer(u0, u0)), v)*dx 
elif CONV == 'n1':
	a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(grad(u0)*u, v)*dx 
	L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx
elif CONV == 'n2':
	a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + inner(grad(u)*u0, v)*dx 
	L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx
elif CONV == 'n3':
	a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx
	L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx - inner(grad(u0)*u0, v)*dx 
elif CONV == 's1':
	a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + 0.5*inner(grad(u0)*u, v)*dx + 0.5*inner(div(outer(u0, u)), v)*dx 
	L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx
elif CONV == 's2':
	a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + 0.5*inner(grad(u)*u0, v)*dx + 0.5*inner(div(outer(u, u0)), v)*dx 
	L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx
elif CONV == 's3':
	a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + 0.5*inner(grad(u0)*u, v)*dx + 0.5*inner(div(outer(u, u0)), v)*dx 
	L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx
elif CONV == 's4':
	a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx + 0.5*inner(grad(u)*u0, v)*dx + 0.5*inner(div(outer(u0, u)), v)*dx 
	L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx
elif CONV == 's5':
	a1 = 1/tau*inner(u, v)*dx + (1/Re)*inner(grad(u), grad(v))*dx
	L1 = inner(f, v)*dx + 1/tau*inner(u0, v)*dx - inner(grad(p0), v)*dx - 0.5*inner(grad(u0)*u0, v)*dx - 0.5*inner(div(outer(u0, u0)), v)*dx 
else:
	print "Unknown CONV", CONV
	exit(-1)
a2 = -inner(grad(p), grad(q))*dx
L2 = 1.0/tau*div(u12)*q*dx - inner(grad(p0), grad(q))*dx
a3 = 1/tau*inner(u, v)*dx
L3 = 1/tau*inner(u12, v)*dx - inner(grad(p1-p0), v)*dx

V_psi = FunctionSpace(mesh, "Lagrange", 2)
psi = TrialFunction(V_psi)
v_psi = TestFunction(V_psi)
psi1 = Function(V_psi)
a_psi = inner(grad(psi), grad(v_psi))*dx
L_psi = inner(-rot(u1), v_psi)*dx
bcs_psi = [DirichletBC(V_psi, Constant(0), domain_all_walls)]

fvelo = File("result_%s_%d_%f_%d/velocity.pvd" % (CONV, N, DT, REYNOLDS))
fpress = File("result_%s_%d_%f_%d/pressure.pvd" % (CONV, N, DT, REYNOLDS))
fpsi = File("result_%s_%d_%f_%d/psi.pvd" % (CONV, N, DT, REYNOLDS))
outpsi = 0
if (MPI.process_number() == 0):
	outpsi = open("result_%s_%d_%f_%d/psi.txt" % (CONV, N, DT, REYNOLDS), "w")
	outpsi.write("t psimax psimin\n")
t = 0
velo_err = 10
while t <= T:
	if (MPI.process_number() == 0):
		print "t = %lf" % t
	solve(a1 == L1, u12, bcs=bcs)
	solve(a2 == L2, p1)
	solve(a3 == L3, u1, bcs=bcs)
	solve(a_psi == L_psi, psi1, bcs=bcs_psi)
	psi_max = psi1.vector().array().max()
	psi_min = psi1.vector().array().min()
	psi_max = MPI.max(psi_max)
	psi_min = MPI.max(psi_min)
	if (MPI.process_number() == 0):
		print "%f %f %f" % (t, psi_max, psi_min)
		outpsi.write("%f %f %f\n" % (t, psi_max, psi_min))
		outpsi.flush()
	fvelo << u1
	fpress << p1
	fpsi << psi1
	u0.assign(u1)
	p0.assign(p1)
	t += DT

#p = plot(u1, mode="color", window_width=1024, window_height=768)
#p.write_png("png/velocity_%s_%d_%f_%d" % (CONV, N, DT, REYNOLDS))
if (MPI.process_number() == 0):
	outpsi.close()