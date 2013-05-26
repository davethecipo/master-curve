#!/usr/bin/python
# -*- coding: latin-1

import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import numpy
from pprint import pprint
from bisect import bisect


def approssima(dato, steps=5000):
  """ Data una serie di punti sperimentali ad una certa temperatura, calcola 
      la spezzata che passa per i punti"""
  tempi = dato['tempi']
  moduli = dato['moduli']
  x_left = numpy.amin(tempi)
  x_right = numpy.amax(tempi)
  y_low = numpy.amin(moduli)
  y_high = numpy.amax(moduli)
  spline =  UnivariateSpline(tempi, moduli, k=1, s=0) # controllare se serve il logaritmo
  risultati = {'x_scale':(x_left,x_right), 'y_scale':(y_low, y_high),'equation':spline}
  return risultati

def draw_points(x,y):
  plt.plot(x,y, 'bo')
  plt.hold('on')
  
def draw_spline(x_left, x_right, spline_equation, draw_steps=100):  
  valori_x_disegno = numpy.linspace(x_left, x_right, draw_steps)
  y_interpolati_disegno = spline_equation(valori_x_disegno)
  plt.plot(valori_x_disegno, y_interpolati_disegno)
  plt.hold('on')    

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
        dati[l]['tempi'].append(numpy.log10(tempo))
        dati[l]['moduli'].append(numpy.log10(modulo))
  
  return dati

def ottieni_informazioni_globali(blocchi_dati):
  maximums = []
  minimums = []
  equations = []
  estremi_modulo = []
  for elem in range(len(blocchi_dati)):
    risultato = approssima(blocchi_dati[elem])
    maximums.append(risultato['x_scale'][1])
    minimums.append(risultato['x_scale'][0])
    equations.append(risultato['equation'])
    estremi_modulo.append(risultato['y_scale'])
  tempo_max = numpy.amax(maximums)
  tempo_min = numpy.amin(minimums)

  return {'equations':equations, 'tempo_max': tempo_max, 'tempo_min':
          tempo_min,'estremi_modulo':estremi_modulo}
  
def calcola_tempi(lista1, lista2, y_min, y_max):
  """ Restituisce l'indice della lista tale per cui i valori dei tempi 
      sono relativi ai punti in cui il modulo assume il valore medio
      rispetto ai valori comuni di modulo.
      
      lista1 è relativo ad una temperatura maggiore di lista2
      
      I valori minimo e massimo sono facilmente determinabili nel caso in cui
      si sappia già che i dati nell'array spline_approssimate sono ordinati
      per temperatura crescente. Per semplicità supponiamo che l'intersezione
      sia non nulla, bisognerebbe controllare intersezioni vuote onde evitare
      interruzioni impreviste del programma      
      """
  modulo_medio = numpy.around((y_max+y_min)/2, 10) # 10 cifre decimali
  indice1 = bisect(lista1, modulo_medio)
  indice2 = bisect(lista2, modulo_medio)
  return indice1, indice2
  
def calcola_a(tempi):
  """ Calcola lo shift factor sapendo il valore del logaritmo dei tempi"""
  counter = 0
  a = [tempi[1]-tempi[0]]
  for i in range(3, len(tempi), 2):
    """ attenzione ai segni: questo ciclo ritorna il valore dello shift factor 
        relativo al passaggio da una temperatura minore ad una maggiore"""
    a_t = tempi[i]-tempi[i-1]+a[counter]
    counter += 1
    a.append(a_t)
  return a
  
  
def ricava_wlf(shift_factors, temperature):
  """ Dati i valori sperimentali di shift factor e temperature, ritorna
      i parametri della WLF"""  
  y = [1/shift_factors[i] for i in range(len(shift_factors))]
  x = []
  for i in range(1, len(temperature)):
    x.append(1/(temperature[i]-temperature[0]))
    
  plt.plot(x,y, 'bo')
  plt.hold('on')
  A = numpy.vstack([x, numpy.ones(len(x))]).T
  m, q = numpy.linalg.lstsq(A, y)[0]
  x_interp = numpy.linspace(x[0],x[-1],50)
  y_interp = m*x_interp+q
  plt.plot(x_interp,y_interp)
  plt.hold('on')
  return m, q

FILE = '/home/davide/poli/2/2-semestre/polimeri/A/Cedevolezza-PS-es5.txt'
CALC_STEPS = 5000

punti_sperimentali = sorted(open_from_csv(FILE), key=lambda t: t['temperatura']) # mette in ordine crescente di temperatura
risultati = ottieni_informazioni_globali(punti_sperimentali)
griglia_x = numpy.linspace(risultati['tempo_min'], risultati['tempo_max'], CALC_STEPS)
spline_approssimate = [ elem(griglia_x) for elem in risultati['equations'] ]
elenco_temperature = [punti_sperimentali[elem]['temperatura'] for elem in range(len(punti_sperimentali)) ]



tempi = []
for i in range(1,len(punti_sperimentali)):
  a, b = calcola_tempi(spline_approssimate[i], spline_approssimate[i-1], risultati['estremi_modulo'][i][0], risultati['estremi_modulo'][i-1][1])
  tempi.append(griglia_x[a])
  tempi.append(griglia_x[b])
  
shift_factor = calcola_a(tempi) # questi sono ricavati sperimentalmente

m, q = ricava_wlf(shift_factor, elenco_temperature)
a = -1/q
b = -a*m
  
for i in range(1,len(elenco_temperature)):
  delta_T = elenco_temperature[i]-elenco_temperature[0]
  log_a = -(a*delta_T)/(b+delta_T)
  # print("teorico da parametri WLF: " + str(log_a))
  nuovi_x = [elem + shift_factor[i-1] for elem in punti_sperimentali[i]['tempi']]
  plt.plot(nuovi_x, punti_sperimentali[i]['moduli'], 'bo')
  plt.hold('on')
"""punti_teorici_x = [punti_sperimentali]
punti_y = numpy.concatenate([punti_sperimentali[i]['tempi'] for i in range(len(punti_sperimentali))])
for elem in range(len(punti_sperimentali)):
"""  


# grafico sperimentale

for i in range(len(punti_sperimentali)):
  for elem in range(len(punti_sperimentali[i]['tempi'])):
    punto_x = punti_sperimentali[i]['tempi'][elem]
    punto_y = punti_sperimentali[i]['moduli'][elem]
    draw_points(punto_x, punto_y)
  eq = risultati['equations'][i]
  xmin = punti_sperimentali[i]['tempi'][0]
  xmax = punti_sperimentali[i]['tempi'][-1]
  draw_spline(xmin, xmax, eq)
    
#for point in range(len()):
  














plt.show()
