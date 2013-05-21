#!/usr/bin/python
# -*- coding: latin-1


import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, splrep, splev
import numpy
from pprint import pprint


    
class experimentalDataBlock():
  """questa classe rappresenta un insieme di dati alla stessa temperatura"""
  def __init__(self, in_temp, in_times, in_modules):
    """in_times e in_modules sono degli array contenenti i valori di tempi e moduli"""
    self.temp = in_temp
    self.times = numpy.array(in_times)
    self.modules = numpy.array(in_modules)
    
  
  
  def logScale(self, in_numpy_array):
    """restituisci i valori in entrata in scala logaritmica"""
    out_array = []
    for i in range(len(in_numpy_array)):
      out_array.append(numpy.log(in_numpy_array[i]))
      
    return out_array
  
  def findSpline(self):
    """trova la spline che passa per tutti i punti sperimentali"""
    s = UnivariateSpline(self.logScale(self.times), self.logScale(self.modules), k=1, s=0)
    return s

  def borders(self, in_numpy_array):
    """trova i valori minimo e massimo della grandezza"""
    min = numpy.amin(in_numpy_array)
    max = numpy.amax(in_numpy_array)
    return min, max
    
    
  def __str__(self):
    print("temperatura della serie di dati: " + str(self.temp))
    print("tempo [s]  ----> modulo ")
    for i in range(len(self.modules)):
      print(str(self.times[i]) + "\t\t" + str(self.modules[i]))
    # resistuisce errore, il problema potrebbe essere l'ultimo elemento della lista vuoto?!


  def drawExperimentalPoints(self):
    plt.plot(self.times, self.modules, 'bo')
    plt.hold('on')
  
  
  def drawSpline(self, steps=10000):
    x_left, x_right = self.borders(self.times)
    print(x_left)
    print(x_right)
    approx_function = self.findSpline()
    x_points = numpy.linspace(x_left, x_right, steps)
    y_interpolated_points = approx_function(x_points)
    plt.plot(x_points, y_interpolated_points)
    plt.hold('on')
    

def open_from_csv(in_file):
  raw = open(in_file,'r', encoding='latin-1').readlines()
  raw = [elem.replace(',','.').replace('\r\n','') for elem in raw]
  raw_temps = raw[0].split('\t')[1:]  
  raw_temps = [elem.lstrip('J(').rstrip('°C)(Pa-1)') for elem in raw_temps]
  last_raw_temp = [raw_temps[-1].rstrip('°C)(Pa-1)\n')]  
  raw_temps = raw_temps[:-1] + last_raw_temp
  raw_times = [elem.split('\t')[0] for elem in raw[1:]]
  
  datamatrix = [[float(temperatura),[],[]] for temperatura in raw_temps]
  
  for i in range(len(raw_times)):
    tempo = float(raw_times[i])
    for l in range(len(raw_temps)):
      dato = raw[i+1].split('\t')[l+1].replace('\n','')
      if dato is not "":
        modulo = float(dato)
        datamatrix[l][1].append(tempo)
        datamatrix[l][2].append(modulo)
  
 # return datamatrix
  blocchi_di_dati = []
  for l in range(len(datamatrix)):
    tempi = datamatrix[l][1]
    moduli = datamatrix[l][2]
    temperatura = datamatrix[l][0]
    dato = experimentalDataBlock(temperatura, tempi, moduli)
    blocchi_di_dati.append(dato)
  return blocchi_di_dati

temps = open_from_csv('/home/davide/poli/2/2-semestre/polimeri/A/Cedevolezza-PS-es5.txt')
#boh = temps[0]

#boh.drawExperimentalPoints()
# è decisamente cannato boh.drawSpline()

buh = temps[1]
buh.drawExperimentalPoints()

for elem in range(len(temps)):
  temps[elem].drawExperimentalPoints()


plt.xscale('log')
plt.yscale('log')
plt.show()