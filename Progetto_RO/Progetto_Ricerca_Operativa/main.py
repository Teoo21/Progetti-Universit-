from instances import *
from environments import *
from solvers import *
from solutions import *

instance_name = 'dummy_problem'

inst = Instance(instance_name)
env = Environment(inst)
solver = solver_307402_322061_323645(env)

X, Y = solver.solve()

sol = Solution(X, Y)
sol.write(instance_name)