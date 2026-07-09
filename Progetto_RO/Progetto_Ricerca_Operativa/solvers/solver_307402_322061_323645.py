from .abstract_solver import AbstractSolver
import numpy as np
from gurobipy import Model, GRB, quicksum
import os

class solver_307402_322061_323645(AbstractSolver):
    def __init__(self, env):
        super().__init__(env)
        self.name = 'solver_307402_322061_323645'
    
    def solve(self):
        super().solve()

        instance = self.env.inst
        service_matrix = instance.service  
        distances = instance.distances  
        weights = instance.weights

        num_magazzini = service_matrix.shape[0]
        num_supermarkets = service_matrix.shape[1]
        num_posizioni = num_magazzini + 1  

        model = Model("PolliTech")

        x = model.addVars(num_magazzini, vtype=GRB.BINARY, name="x")
        y = model.addVars(num_posizioni, num_posizioni, vtype=GRB.BINARY, name="y")
        u = model.addVars(num_posizioni, vtype=GRB.CONTINUOUS, lb=0, ub=num_posizioni, name="u")
        z_s = model.addVars(num_supermarkets, vtype=GRB.BINARY, name="z_s")

        penalità = quicksum(weights['missed_supermarket'] * z_s[s] for s in range(num_supermarkets))
        
        for s in range(num_supermarkets):
            model.addConstr(quicksum(service_matrix[m][s] * x[m] for m in range(num_magazzini)) >= 1 - z_s[s])

        model.setObjective(
            weights['construction'] * quicksum(x[m] for m in range(num_magazzini)) +
            penalità +
            weights['travel'] * quicksum(distances[i][j] * y[i, j] for i in range(num_posizioni) for j in range(num_posizioni)),
            GRB.MINIMIZE)

        for i in range(1, num_posizioni):
            model.addConstr(quicksum(y[i, v] for v in range(num_posizioni) if v != i) == x[i - 1])

        for v in range(1, num_posizioni):
            model.addConstr(quicksum(y[u, v] for u in range(num_posizioni) if u != v) == x[v - 1])

        model.addConstr(quicksum(y[0, j] for j in range(1, num_posizioni)) == 1)
        model.addConstr(quicksum(y[j, 0] for j in range(1, num_posizioni)) == 1)

        model.addConstr(u[0] == 0.0)  
        for u_i in range(1, num_posizioni):
            model.addConstr(u[u_i] >= 2)
            model.addConstr(u[u_i] <= num_posizioni)

        for u_i in range(1, num_posizioni):
            for v_i in range(1, num_posizioni):
                if u_i != v_i:
                    model.addConstr(u[u_i] - u[v_i] + num_posizioni * y[u_i, v_i] <= num_posizioni - 1)

        model.setParam('OutputFlag', 0)
        model.optimize()

        X = [int(x[i].X >= 0.5) for i in range(num_magazzini)]
        Y = [[int(y[i, j].X >= 0.5) for j in range(num_posizioni)] for i in range(num_posizioni)]

        return X, Y