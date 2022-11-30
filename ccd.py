#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Robótica Computacional - 
# Grado en Ingeniería Informática (Cuarto)
# Práctica: Resolución de la cinemática inversa mediante CCD
#           (Cyclic Coordinate Descent).

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import colorsys as cs

# ******************************************************************************
# Declaración de funciones

def muestra_origenes(O, final = 0):
  # Muestra los orígenes de coordenadas para cada articulación
  print('Origenes de coordenadas:')
  for i in range(len(O)):
    print('(O' + str(i) + ')0\t= ' + str([round(j, 3) for j in O[i]]))
  if final:
    print('E.Final = ' + str([round(j, 3) for j in final]))

def muestra_robot(O, obj, data):
  # Muestra el robot graficamente
  plt.figure(1)
  plt.xlim(-L, L)
  plt.ylim(-L, L)
  T = [np.array(o).T.tolist() for o in O]
  for i in range(len(T)):
    plt.plot(T[i][0], T[i][1], '-o', color = cs.hsv_to_rgb(i / float(len(T)), 1, 1))
    for j in range(len(T[i])):
      if j < len(data) - 1:
        if data[j]["type"] == 'p':
          plt.plot(T[i][0][j], T[i][1][j], '-s', color = cs.hsv_to_rgb(i / float(len(T)), 1, 1), ms=7)
  plt.plot(obj[0], obj[1], '*')
  plt.show()
  raw_input()
  plt.clf()

# def muestra_robot(O,obj):
# # Muestra el robot graficamente
# plt.figure(1)
# plt.xlim(-L,L)
# plt.ylim(-L,L)
# T = [np.array(o).T.tolist() for o in O]
# for i in range(len(T)):
#   plt.plot(T[i][0], T[i][1], '-o', color=cs.hsv_to_rgb(i/float(len(T)),1,1))
# plt.plot(obj[0], obj[1], '*')
# plt.show()
# raw_input()
# plt.clf()

def matriz_T(d, th, a, al):
  # Calcula la matriz T (ángulos de entrada en grados)
  
  return [[cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th)]
         ,[sin(th),  cos(th)*cos(al), -sin(al)*cos(th), a*sin(th)]
         ,[      0,          sin(al),          cos(al),         d]
         ,[      0,                0,                0,         1]
         ]

def cin_dir(th,a):
  #Sea 'th' el vector de thetas
  #Sea 'a'  el vector de longitudes
  T = np.identity(4)
  o = [[0, 0]]
  for i in range(len(th)):
    T = np.dot(T,matriz_T(0, th[i], a[i], 0))
    tmp = np.dot(T, [0, 0, 0, 1])
    o.append([tmp[0], tmp[1]])
  return o

# ******************************************************************************
# Cálculo de la cinemática inversa de forma iterativa por el método CCD

# valores articulares arbitrarios para la cinemática directa inicial
#th = [0., 0., 0.]
#a = [5., 5., 5.]
#L = sum(a) # variable para representación gráfica
EPSILON = .01

plt.ion() # modo interactivo

# introducción del punto para la cinemática inversa
if len(sys.argv) != 6:
  sys.exit("python " + sys.argv[0] + "input.txt límite_prismático limite_rotación punto_x punto_y")


data = []
with open(sys.argv[1]) as file:
  for line in file:
    split = line.split((' '))
    data.append({
      "th": float(split[1]),
      "a": float(split[0]),
      "type": split[2][0]
    })

objetivo = [float(i) for i in sys.argv[4:]]

prismatic_limit = int(sys.argv[2])
rotation_limit = float(sys.argv[3]) * np.pi / 180

L = sum([art["a"] for art in data])

O = range(len(data) + 1) # Reservamos estructura en memoria
O[0] = cin_dir([art["th"] for art in data], [art["a"] for art in data]) # Posicion inicial
print("- Posicion inicial:")
muestra_origenes(O[0])

dist = float("inf")
prev = 0.
iteracion = 1
while (dist > EPSILON and abs(prev - dist) > EPSILON / 100.):
  prev = dist
  # Para cada combinación de articulaciones:
  for i in range(len(data)):
    # Cálculo de la cinemática inversa:
    # Articulaciones de rotación
    if data[i]["type"] == 'r':
      v1  = [O[i][-1][0] - O[i][-i - 2][0], O[i][-1][1] - O[i][-i - 2][1]]
      v2 = [objetivo[0] - O[i][-i - 2][0], objetivo[1] - O[i][-i - 2][1]]
      
      alfa2 = atan2(v1[1], v1[0])
      alfa1 = atan2(v2[1], v2[0])

      inc_th = alfa1 - alfa2

      # correción de ángulo
      if inc_th > np.pi:                 # Rotación > 180º
        inc_th = 0 - (2*np.pi - inc_th)  # 360 - 270 = 90 -> -90
      elif inc_th < -np.pi:              # Rotación < -180º
        inc_th = 2*np.pi + inc_th        # 360 + (-270) = 90
      
      # límite de rotación
      if data[len(data) - i - 1]["th"] + inc_th < rotation_limit:
        data[len(data) - i - 1]["th"] = data[len(data) - i - 1]["th"] + inc_th
      else:
        data[len(data) - i - 1]["th"] = rotation_limit

    # Articulaciones prismáticas
    if data[i]["type"] == 'p':     
      w = sum([art["th"] for art in data[0:i]])
      
      resta = [objetivo[0] - O[i][-1][0], objetivo[1] - O[i][-1][1]]
      
      # límite prismática
      if data[i]["a"] + cos(w) * resta[0] + sin(w) * resta[1] < prismatic_limit:
        data[i]["a"] += cos(w) * resta[0] + sin(w) * resta[1]
      else:
        data[i]["a"] = prismatic_limit
      L = sum([art["a"] for art in data])
    O[i + 1] = cin_dir([art["th"] for art in data], [art["a"] for art in data])

  dist = np.linalg.norm(np.subtract(objetivo, O[-1][-1]))
  print("\n- Iteracion " + str(iteracion) + ':')
  muestra_origenes(O[-1])
  muestra_robot(O, objetivo, data)
  print("Distancia al objetivo = " + str(round(dist, 5)))
  iteracion += 1
  O[0] = O[-1]

if dist <= EPSILON:
  print("\n" + str(iteracion) + " iteraciones para converger.")
else:
  print("\nNo hay convergencia tras " + str(iteracion) + " iteraciones.")

print("- Umbral de convergencia epsilon: " + str(EPSILON))
print("- Distancia al objetivo:          " + str(round(dist, 5)))
print("- Valores finales de las articulaciones:")
for i in range(len(data)):
  print("  theta" + str(i + 1) + " = " + str(round(data[i]["th"], 3)))
for i in range(len(data)):
  print("  L" + str(i + 1) + "     = " + str(round(data[i]["a"], 3)))
