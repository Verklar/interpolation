import math
import numpy as np
import scipy.interpolate

#функция задана таблицей
x1 = np.array([2.6, 3.3, 4.7, 6.1, 7.5, 8.2, 9.6])
y1 = np.array([2.1874, 2.8637, 3.8161, 3.8524, 3.1905, 2.8409, 2.6137])


def lagrange(a, x, y):
  sum = 0
  for i in range (0, len(y)):
    new = y[i]
    for j in range (0, len(x)):
      if (i!=j): new*=(a - x[j])/(x[i] - x[j])
    sum += new
  return sum

def f(x):
  return 1.7 * pow(x, 1/3) - math.cos(0.4 - 0.7*x)

print('Проверка, совпадает ли ИМЛ с табличными значениями')
for i in range (0, len(x1)):
  result = lagrange(x1[i], x1, y1)
  print(x1[i], result, result - y1[i] <= 10**-4)


#функция задана таблицей, значения с равноотстающим шагом
x2 = np.array([2.6, 3.6, 4.6, 5.6, 6.6, 7.6, 8.6])
y2 = np.array([2.187, 3.127, 3.776, 3.948, 3.662, 3.136, 2.695])


def Newton_1(a, x, y):
  n = len(x2)
  h = x2[1] - x2[0]
  f = np.zeros((7, 7))#таблица конечных разностей
  for i in range(0, n):
      f[0][i] = y2[i]
  for k in range(1, n):
      for i in range(n - k):
          f[k][i] = (f[k - 1][i + 1] - f[k - 1][i])
  np.set_printoptions(precision=3, suppress=True, formatter={'all': lambda x: f'{x:0.3f}'})#вывод с тремя знаками после запятой
  f = f.transpose()

  N1 = f[0][0]
  for i in range (1, n):
    mult = f[0][i]/(math.factorial(i)*pow(h, i))
    for j in range (i):
      mult*=(a-x[j])
    N1+=mult
  return(N1)



def Newton_2(a, x, y):
  n = len(x2)
  h = x2[1] - x2[0]
  f = np.zeros((7, 7))#таблица конечных разностей
  for i in range(0, n):
      f[0][i] = y2[i]
  for k in range(1, n):
      for i in range(n - k):
          f[k][i] = (f[k - 1][i + 1] - f[k - 1][i])
  np.set_printoptions(precision=3, suppress=True, formatter={'all': lambda x: f'{x:0.3f}'})#вывод с тремя знаками после запятой
  f = f.transpose()

  N1 = f[n-1][0]
  for i in range (1, n):
    mult = f[n-i-1][i]/(math.factorial(i)*pow(h, i))
    for j in range (i):
      mult*=(a-x[n - 1 -j])
    N1+=mult
  return(N1)

print('Проверим соответствие значений интерполяций и функции для точки 4.0')

a = 4.0
res1 = lagrange(a, x1, y1)
print('ИМЛ для числа 4.0 -- ', res1, 'разница при вычислении по функции и по ИМЛ -- ', abs(f(a)-res1))

res2 = Newton_1(a, x2, y2)
print('Первый ИМН для числа 4.0 -- ', res2, 'разница при вычислении по функции и по первому ИМН -- ', abs(f(a)-res2))


res3 = Newton_2(a, x2, y2)
print('Второй ИМН для числа 4.0 -- ', res3, 'разница при вычислении по функции и по второму ИМН -- ', abs(f(a)-res3))

print('Проверим соответствие значений стандартных интерполяций и функции для точки 4.0')

g = lagrange(a, x2, y2)
print('scipy ИМЛ для числа 4.0 -- ', g, 'разница при вычислении по функции и по ИМЛ -- ', abs(f(a)-g))

cs = scipy.interpolate.CubicSpline(x2, y2)
csres = cs(a)
print('scipy cubicspline для числа 4.0 -- ', csres, 'разница при вычислении по функции и cubicspline -- ', abs(f(a)-csres))

pc = scipy.interpolate.PchipInterpolator(x2, y2)
pcres = pc(a)
print('scipy PchipInterpolator для числа 4.0 -- ', pcres, 'разница при вычислении по функции и PchipInterpolator -- ', abs(f(a)-pcres))

cubic = scipy.interpolate.Akima1DInterpolator(x2, y2)
akimares = cubic(a)
print('scipy Akima1DInterpolator для числа 4.0 -- ', akimares, 'разница при вычислении по функции и Akima1DInterpolator -- ', abs(f(a)-akimares))

## лучшие результаты у первого и второго ИМН