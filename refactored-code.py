#!/usr/bin/python
# -*- coding: latin-1


import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import numpy
from pprint import pprint
from bisect import bisect


def approssima(dato, steps=5000):
  x_left = numpy.amin(dato['tempi'])
  x_right = numpy.amax(dato['tempi'])
  spline =  UnivariateSpline(dato['tempi'], dato['moduli'], k=1, s=0)
  y_min = numpy.amin(dato['moduli']) 
  y_max = numpy.amax(dato['moduli'])
  risultati = {'x_scale':(x_left,x_right), 'y_scale':(y_min,y_max),'equation':spline}
  return risultati



def trova_vicino(lista_target, lista_riferimento):
  n=0
  indexes = []
  for i in range(len(lista_riferimento)):
    if lista_target[n] >= lista_riferimento[0] and lista_target[n] <= lista_riferimento[-1]: 
      nearest_found = False
      while nearest_found == False:
        diff = lista_riferimento[i] -  lista_target[n]
        next_diff = lista_riferimento[i] - lista_target[n+1]
        if next_diff > diff:
          nearest_found = True
          indexes.append(n)
        else:
          n += 1
          nearest_found = False
          
  return indexes


""" da riscrivere....  
def draw(x_left, x_right, draw_steps=500, spline_equation):  
  valori_x_disegno = numpy.linspace(x_left, x_right, draw_steps)
  y_interpolati_disegno = spline_equation(valori_x_disegno)
  
  plt.plot(dato['tempi'], dato['moduli'], 'bo')
  plt.hold('on')
  plt.plot(valori_x_disegno, y_interpolati_disegno)
  plt.hold('on')
"""    

def trova_intersezione(intervallo_riferimento, intervallo):
  min, max = 0,0
  if intervallo_riferimento[0] >= intervallo[0]:
    if intervallo_riferimento[1] <= intervallo[1]:
      min, max = intervallo_riferimento[0], intervallo_riferimento[1]
    else:
      min, max = intervallo_riferimento[0], intervallo[1]
      
  else:
    if intervallo_riferimento[1] <= intervallo[1]:
      min, max = intervallo[0], intervallo_riferimento[1]
    else:
      min, max = intervallo[0], intervallo[1]
      
  return min, max
  
    

def open_from_csv(in_file):
  raw = open(in_file,'r', encoding='latin-1').readlines()
  raw = [elem.replace(',','.').replace('\r\n','') for elem in raw]
  raw_temps = raw[0].split('\t')[1:]  
  raw_temps = [elem.lstrip('J(').rstrip('°C)(Pa-1)') for elem in raw_temps]
  last_raw_temp = [raw_temps[-1].rstrip('°C)(Pa-1)\n')]
  
  
  raw_temps = raw_temps[:-1] + last_raw_temp
  raw_times = [elem.split('\t')[0] for elem in raw[1:]]
  
  dati = [{'temperatura':float(temperatura),'moduli':[],'tempi':[]} for temperatura in raw_temps]
  
  for i in range(len(raw_times)):
    tempo = float(raw_times[i])
    for l in range(len(raw_temps)):
      dato = raw[i+1].split('\t')[l+1].replace('\n','')
      if dato is not "":
        modulo = float(dato)
        dati[l]['tempi'].append(tempo)
        dati[l]['moduli'].append(modulo)
  
  return dati



def ottieni_informazioni_globali(blocchi_dati):
  maximums = []
  minimums = []
  equations = []
  for elem in range(len(blocchi_dati)):
    risultato = approssima(blocchi_dati[elem])
    maximums.append(risultato['x_scale'][1])
    minimums.append(risultato['x_scale'][0])
    equations.append(risultato['equation'])
  max_assoluto = numpy.amax(maximums)
  min_assoluto = numpy.amin(minimums)
  
  return {'equations':equations, 'max_assoluto': max_assoluto, 'min_assoluto':
          min_assoluto}
  


FILE = '/home/davide/poli/2/2-semestre/polimeri/A/Cedevolezza-PS-es5.txt'

punti_sperimentali = sorted(open_from_csv(FILE), key=lambda t: t['temperatura']) # mette in ordine crescente di temperatura
risultati = ottieni_informazioni_globali(punti_sperimentali)











riferimento = punti_sperimentali[0] 
risultati_riferimento = approssima(riferimento)

steps = 5000




def trova_tempi(indici, lista):
  valori_x = []
  for i in range(len(indici)):
    valori_x.append(lista[indici[i]])
  return valori_x
    
 

def regressione(riferimento, target):
  risultati = approssima(target)
  min_y, max_y = trova_intersezione(risultati['y_scale'], risultati_riferimento['y_scale']) #risultati_riferimento è una variabile globale, BAD
  low_index = risultati['spline'].index(min_y) if min_y in risultati['spline'] else None
  high_index = risultati['spline'].index(max_y) if max_y in risultati['spline'] else None
  
   
 helper_calcolo = numpy.linspace(min_y, max_y, steps)
  indici_target = trova_vicino(risultati, helper_calcolo)
  indici_riferimento = trova_vicino(dati_riferimento, helper_calcolo) # indice target
  tempi_target = trova_tempi(indici_target, dati_target)
  tempi_riferimento = trova_tempi(indici_riferimento, dati_riferimento)
  a_t = [1/(numpy.log(e[1]/e[0])) for e in zip(tempi_riferimento, tempi_target)]
  x = [1/(riferimento['temperatura']-target['temperatura']) for elem in a_t]
  return x, a_t
  
  
punti_x, punti_y = [],[]

for i in range(1,len(asd)):
  x, y = regressione(riferimento, asd[i])
  punti_x.append(x)
  punti_y.append(y)
  
  
# regressione lineare
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  