#!/usr/bin/python
# -*- coding: UTF-8

import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy import special
import numpy
from pprint import pprint
from bisect import bisect

FILE = '/home/davide/poli/2/2-semestre/polimeri/A/Cedevolezza-PS-es5.txt'
CALC_STEPS = 5000
DRAW_STEPS = 100
REQUIRED_TIME = 315360000 # tempo in secondi al quale vogliamo un certo modulo
OUTPUT_NAME = 'results.txt'
ENTALPIA = 51500*4.184 # J/mol
R = 8.314 # J/mol K
def approssima(dato, steps=CALC_STEPS):
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
  
def draw_spline(x_left, x_right, spline_equation, in_label, draw_steps=DRAW_STEPS):  
  valori_x_disegno = numpy.linspace(x_left, x_right, draw_steps)
  y_interpolati_disegno = spline_equation(valori_x_disegno)
  plt.plot(valori_x_disegno, y_interpolati_disegno, label=in_label)
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
        relativo al passaggio da una temperatura maggiore ad una minore"""
    a_t = tempi[i]-tempi[i-1]+a[counter]
    counter += 1
    a.append(a_t)
  return a

def arrhenius(ent, T,T_rif):
  return numpy.exp(ent/R*(1/T-1/T_rif))

  
  
def ricava_wlf(shift_factors, temperature):
  """ Dati i valori sperimentali di shift factor e temperature, ritorna
      i parametri della WLF"""  
  y = [1/shift_factors[i] for i in range(len(shift_factors))]
  x = []
  for i in range(1, len(temperature)):
    x.append(1/(temperature[i]-temperature[0]))
  A = numpy.vstack([x, numpy.ones(len(x))]).T
  m, q = numpy.linalg.lstsq(A, y)[0]
  x_interp = numpy.linspace(x[0],x[-1],50)
  y_interp = m*x_interp+q
  
  return m, q, [x, y, x_interp, y_interp]

punti_sperimentali = sorted(open_from_csv(FILE), key=lambda t: t['temperatura']) # mette in ordine crescente di temperatura
risultati = ottieni_informazioni_globali(punti_sperimentali)
griglia_x = numpy.linspace(risultati['tempo_min'], risultati['tempo_max'], CALC_STEPS)
spline_approssimate = [ elem(griglia_x) for elem in risultati['equations'] ] # tutte le equazioni della spline
elenco_temperature = [punti_sperimentali[elem]['temperatura'] for elem in range(len(punti_sperimentali)) ] # elenco di temperature
tempi = []
for i in range(1,len(punti_sperimentali)):
  a, b = calcola_tempi(spline_approssimate[i], spline_approssimate[i-1], risultati['estremi_modulo'][i][0], risultati['estremi_modulo'][i-1][1])
  tempi.append(griglia_x[a])
  tempi.append(griglia_x[b])
  
shift_factor = calcola_a(tempi) # questi sono ricavati sperimentalmente
m, q, dati_wlf = ricava_wlf(shift_factor, elenco_temperature)
a = -1/q
b = -a*m


tempi_totali = [punti_sperimentali[0]['tempi']]
for i in range(1,len(elenco_temperature)):
  tempi_totali.append([elem + shift_factor[i-1] for elem in punti_sperimentali[i]['tempi']])

"""penso che si possa evitare di usare sort, ci sto pensando"""
tempi_totali = numpy.sort(numpy.concatenate(tempi_totali))
moduli_totali = numpy.sort(numpy.concatenate([punti_sperimentali[elem]['moduli'] for elem in range(len(punti_sperimentali))]))

eq_master_sperimentale = UnivariateSpline(tempi_totali, moduli_totali, k=1, s=0)
x_totali = numpy.linspace(tempi_totali[0], tempi_totali[-1], CALC_STEPS)
y_interpolati_sper = eq_master_sperimentale(x_totali)
pprint(numpy.log10(arrhenius(ENTALPIA, 25+273, 97+273)))
pprint( numpy.log10(REQUIRED_TIME))

required_time = numpy.log10(REQUIRED_TIME)-numpy.log10(arrhenius(ENTALPIA, 25+273, 97+273)) # 10 anni corretti da shift factor
print(required_time)
indice_modulo = bisect(x_totali, required_time)
pprint(indice_modulo)
modulo_shiftato = special.exp10(y_interpolati_sper[indice_modulo])
pprint(y_interpolati_sper[indice_modulo])

nu = 0.3
f = 0.01 # m
L = 0.3 # m
P = 2940 # Pa



spessore = special.cbrt((3*L**4*P*modulo_shiftato)/(4*f*(1+nu)))
pprint(spessore) # m

# Inizia il codice per i grafici




f1 = plt.figure() 
plt.title("Master curve")
plt.xlabel(r"$log(t) \,\, [s]$")
plt.ylabel(r"$log(J) \,\, [Pa^{-1}]$")
plt.plot(x_totali, y_interpolati_sper, 'bo')
plt.hold('on')

f2 = plt.figure()


#
# grafico sperimentale
#

plt.subplot(221)
plt.title("Dati sperimentali")
plt.xlabel(r"$log(t) \,\, [s]$")
plt.ylabel(r"$log(J) \,\, [Pa^{-1}]$")
for i in range(len(punti_sperimentali)):
  for elem in range(len(punti_sperimentali[i]['tempi'])):
    punto_x = punti_sperimentali[i]['tempi'][elem]
    punto_y = punti_sperimentali[i]['moduli'][elem]
    draw_points(punto_x, punto_y)
  eq = risultati['equations'][i]
  xmin = punti_sperimentali[i]['tempi'][0]
  xmax = punti_sperimentali[i]['tempi'][-1]
  etichetta = str(punti_sperimentali[i]['temperatura']) + " °C"
  draw_spline(xmin, xmax, eq, etichetta)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# WLF
plt.subplot(222)
plt.title("Regressione lineare (WLF)")
plt.xlabel(r"$(T-T_{0})^{-1} \,\, [\degree C^{-1}]$")
plt.ylabel(r"$(log(a_{T}^{T_0}))^{-1} \,\, [s^{-1}]$")
plt.plot(dati_wlf[0], dati_wlf[1], 'bo')
plt.hold('on')
plt.plot(dati_wlf[2],dati_wlf[3])
plt.hold('on')



  
#
# master curve 
#
for i in range(1,len(elenco_temperature)):
  nuovi_x = [elem + shift_factor[i-1] for elem in punti_sperimentali[i]['tempi']]
  
  plt.subplot(223)
  plt.title("Master curve (shift factors sperimentali)")
  plt.xlabel(r"$log(t) \,\, [s]$")
  plt.ylabel(r"$log(J) \,\, [Pa^{-1}]$")
  plt.plot(nuovi_x, punti_sperimentali[i]['moduli'], 'bo')
  plt.hold('on')

  delta_T = elenco_temperature[i]-elenco_temperature[0]
  log_a = -(a*delta_T)/(b+delta_T)
  x_teorici = [elem + log_a for elem in punti_sperimentali[i]['tempi']]

  plt.subplot(224)
  plt.title("Master curve (shift factors da WLF)")
  plt.xlabel(r"$log(t) \,\, [s]$")
  plt.ylabel(r"$log(J) \,\, [Pa^{-1}]$")
  plt.plot(x_teorici, punti_sperimentali[i]['moduli'])
  plt.hold('on')

plt.tight_layout()
plt.show()


output_file = open(OUTPUT_NAME,'w')
output_file.write('Parametri WLF: A=' + str(a) + ', B=' + str(b) + "\n")
output_file.write('Modulo [1/Pa] dopo 10 anni: ' + str(modulo_shiftato) + "\n" )
output_file.write("shift factors sperimentali: \n")
for elem in range(len(shift_factor)):
  output_file.write(str(elenco_temperature[elem+1]) + " °C => " + str(shift_factor[elem])+ "\n")
output_file.close()
