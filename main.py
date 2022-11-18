from numpy import sqrt, arctan, random
import numpy
from math import pi
from matplotlib.pyplot import plot, show, title, xlabel, ylabel

class Atmosphere():
  screen50mm = sqrt(24**2 + 36**2)
  nbPoints = 21
  diapo, screen, screenTarget, diapoCorrected, screenFinal = [], [], [], [], []
  __sizeX = 36
  __sizeY = 24
  __stop = 10**-11

  def __init__(self, g:float, alpha:float, beta:float, focal:float, error:float, nbPixels:int):
    self.g = g
    self.alpha = alpha
    self.beta = beta
    self.nbPixels = nbPixels
    self.gamma2 = 1 + alpha**2
    self.gamma = sqrt(self.gamma2)
    self.focal = focal
    self.error = 1 / error if error != 0 else 0
    self.gM1 = 0
    self.alphaM1 = 0
    if not self.diapo: self.createDiapo()
    if not self.screen: self.createScreen()
    self.computeParameters()
    if not self.screenTarget: self.createScreenTarget()
    if not self.diapoCorrected: self.createDiapoCorrected()
    if not self.screenFinal: self.createScreen(True)
    self.computeGAlpha()
    self.computeBeta()

  def __str__(self):
    return f"Atm g: {self.g}, alpha: {round(arctan(self.alpha) * 180 / pi)}°, beta: {round(arctan(self.beta) * 180 / pi)}°"

  def prettyPrint(self):
    print(f"Atm g: {self.g}, alpha: {round(arctan(self.alpha) * 180 / pi)}°, beta: {round(arctan(self.beta) * 180 / pi)}°")
    print(f"precision : {1/self.error}, nb pixels in one line : {self.nbPixels}")
    print(f"precision g : {self.computePrecision(self.g, self.gGd)}, precision alpha : {self.computePrecision(self.alpha, self.alphaGd)} precision beta : {self.computePrecision(self.beta, self.betaGd)}")
    arrayPixels = [0] * 11
    for i in range(len(self.screenFinal)):
      for j in range(len(self.screenFinal)):
        error = sqrt((self.screenFinal[i][j][0] - (self.g * self.diapo[i][j][0]))**2 + (self.screenFinal[i][j][1] - (self.g * self.diapo[i][j][1]))**2)
        error *= self.nbPixels / self.g / self.sizeX * self.focal
        if error > 1: arrayPixels[10] += 1
        elif error > 0.9: arrayPixels[9] += 1
        elif error > 0.8: arrayPixels[8] += 1
        elif error > 0.7: arrayPixels[7] += 1
        elif error > 0.6: arrayPixels[6] += 1
        elif error > 0.5: arrayPixels[5] += 1
        elif error > 0.4: arrayPixels[4] += 1
        elif error > 0.3: arrayPixels[3] += 1
        elif error > 0.2: arrayPixels[2] += 1
        elif error > 0.1: arrayPixels[1] += 1
        else: arrayPixels[0] += 1
    lower, nbPixels = 0.1, len(self.diapo)**2
    for value in arrayPixels:
      print(f"nombres de pixel de précision supérieure à {lower} : {int(value * 100 / nbPixels)}%")
      lower += 0.1


  def computePrecision(self, value, approx):
    if value - approx != 0:
      if value != 0: return int(abs(value / (value - approx)))
      return int(abs(approx / (value - approx)))
    return "infinity"

  def drawScreen (self, type:str):
      xs, ys = [], []
      if type == "diapo": plots = self.diapo
      elif type == "screen": plots = self.screen
      elif type == "target": plots = self.screenTarget
      elif type == "corrected": plots = self.diapoCorrected
      elif type == "final": plots = self.screenFinal
      else:
        print(f"type {type} does not exist")
        return
      for index in range(self.nbPoints):
        xs += [p[0] for p in plots[index]]
        ys += [p[1] for p in plots[index]]
      plot(xs, ys, 'o', color='black')
      title(type)
      xlabel("X Axis")
      ylabel("Y Axis")
      show()

  @property
  def sizeX(cls): return sqrt(cls.__sizeX**2 + cls.__sizeY**2)

  @property
  def sizeY(cls): return sqrt(cls.__sizeX**2 + cls.__sizeY**2) 

  def createDiapo(self):
    if not self.diapo:
      startX = - self.sizeX / (2 * self.focal)
      stepX = -startX / ((self.nbPoints - 1) / 2)
      startY = - self.sizeY / (2 * self.focal)
      stepY = -startY / ((self.nbPoints - 1) / 2)
      valueY = startY
    for indexY in range(self.nbPoints):
      newColumn, valueX = [], startX
      Atmosphere.diapo.append(newColumn)
      for indexX in range(self.nbPoints):
        newColumn.append((valueX, valueY))
        valueX += stepX
      valueY += stepY

  def createScreen (self, corrected=False, error=True):
    cst = self.g/self.gamma
    numY = sqrt(self.beta**2 + self.gamma2)
    if corrected:
      self.screenFinal = []
    else:
      self.screen = []
    diapo = self.diapoCorrected if corrected else self.diapo
    for line in diapo:
      screenLine = []
      for xy in line:
        numX = self.gamma2 * xy[0] + self.alpha*self.beta*xy[1]
        den = 1 + self.alpha * xy[0] + self.beta * xy[1]
        errorX = random.normal(0, self.error) * self.sizeX / self.focal if error else 0
        errorY = random.normal(0, self.error) * self.sizeY / self.focal  if error else 0
        screenX =  cst * numX / den + errorX
        screenY = cst * numY * xy[1] / den + errorY
        screenLine.append((screenX, screenY))
      if corrected:
        self.screenFinal.append(screenLine)
      else:
        self.screen.append(screenLine)


  def createScreenTarget(self):
    startX = - self.gGd * self.sizeX / (2 * self.focal)
    stepX = -startX / ((self.nbPoints - 1) / 2)
    startY = - self.gGd * self.sizeY / (2 * self.focal)
    stepY = -startY / ((self.nbPoints - 1) / 2)
    valueY = startY
    for indexY in range(self.nbPoints):
      newColumn, valueX = [], startX
      Atmosphere.screenTarget.append(newColumn)
      for indexX in range(self.nbPoints):
        newColumn.append((valueX, valueY))
        valueX += stepX
      valueY += stepY

  def createDiapoCorrected(self):
    gamma2 = 1 + self.alphaGd**2
    gamma = sqrt(gamma2)
    tNorm = sqrt(gamma2 + self.betaGd**2)
    nxx = tNorm
    nxy = - self.alphaGd * self.betaGd
    dc = self.gGd * gamma * tNorm
    dx = -self.alphaGd * tNorm
    dy = - self.betaGd
    nyy = gamma2
    self.diapoCorrected = []
    for line in self.screenTarget:
      screenLine = []
      for xy in line:
        numX = nxx * xy[0] + nxy * xy[1]
        den = dc + dx * xy[0] + dy * xy[1]
        numY = nyy * xy[1]
        screenLine.append((numX / den, numY / den))
      self.diapoCorrected.append(screenLine)

  def computeGAlpha(self):
    A, B, C, D, E = 0, 0, 0, 0, 0
    for i in range(len(self.diapo)):
      for j in range(len(self.diapo)):
        A += self.diapo[i][j][0]**2 * self.screen[i][j][0]**2
        B += self.diapo[i][j][0]**2 * self.screen[i][j][0]
        C += self.diapo[i][j][0]**2
        D += self.diapo[i][j][0] * self.screen[i][j][0]**2
        E += self.diapo[i][j][0] * self.screen[i][j][0]
    self.alphaM1 = (C*D-B*E) / (B**2 - A*C)
    self.gM1 = (1 + self.alphaM1**2)**(-1/2) * (B*D-A*E) / (B**2 - A*C)

  def computeBeta(self):
    A, B, C, D, E, F, G, H = 0, 0, 0, 0, 0, 0, 0, 0
    for i in range(len(self.diapo)):
      for j in range(len(self.diapo)):
        A += self.diapo[i][j][1]**2 * self.screen[i][j][1]**2
        B += self.diapo[i][j][1]**2
        C += self.diapo[i][j][1] * self.screen[i][j][1]**2
        D += self.diapo[i][j][0]*self.diapo[i][j][1] * self.screen[i][j][1]**2
        E += self.screen[i][j][1]**2
        F += self.diapo[i][j][0] * self.screen[i][j][1]**2
        G += self.diapo[i][j][0]**2 * self.screen[i][j][1]**2

    a = A - self.gM1**2 / (1 + self.alphaM1**2) * B
    b = 2*( C + self.alphaM1 * D)
    c = E + 2*self.alphaM1*F + self.alphaM1**2 * G - self.gM1**2 * B
    delta = sqrt(b**2 - 4*a*c)
    self.betaM1 = (-b - delta)/ (2*a)
    if self.betaM1 < 0: beta = (-b + delta)/ (2*a)

  def computeParameters(self):
    paramComputed = self.gradientDescent(self.costFunction, [0,0,0], 1)
    Atmosphere.gGd = paramComputed[0]
    Atmosphere.alphaGd = paramComputed[1]
    Atmosphere.betaGd = paramComputed[2]
    return paramComputed

  def gradientDescent(self, costFunction, array, lambdaVal, cost="unknown", gradient="unknown"):
    if cost == "unknown": cost = costFunction(array)
    if gradient == "unknown": gradient = self.computeGradient(costFunction, array)
    newArray = [array[index] - lambdaVal * gradient[index] for index in range(len(array))]
    newCost = costFunction(newArray)
    if sum([(newArray[index]- array[index])**2 for index in range(len(array))]) < self.__stop:
      return array
    if newCost < cost:
      return self.gradientDescent(costFunction, newArray, lambdaVal=lambdaVal, cost=newCost)
    else:
      return self.gradientDescent(costFunction, array, lambdaVal=lambdaVal/2, cost=cost, gradient=gradient)

  def computeGradient(self, costFunction, array):
    gradient = [0] * len(array)
    for index in range(len(array)):
      arrayPlusDelta = [0] * len(array)
      for subIndex in range(len(array)):
        arrayPlusDelta[subIndex] = array[subIndex] + self.__stop if index == subIndex else array[subIndex]
      gradient[index] = (costFunction(arrayPlusDelta) - costFunction(array)) / self.__stop
    return gradient

  def costFunction (self, array):
    g,alpha,beta= array[0], array[1], array[2]
    gamma2 = 1 + alpha**2
    gamma = sqrt(gamma2)
    cost = 0
    for i in range(len(self.diapo)):
      line = self.diapo[i]
      for j in range(len(line)):
        phiX = self.phiX(line[j], self.screen[i][j], g, alpha, beta, gamma2, gamma)
        phiY = self.phiY(line[j], self.screen[i][j], g, alpha, beta, gamma2, gamma)
        cost += phiX**2 + phiY**2
    return cost

  def phiX(self, xy, xyIm, g, alpha, beta, gamma2, gamma):
    return gamma * (1 + alpha*xy[0] + beta*xy[1]) * xyIm[0] - g * (gamma2 * xy[0] + alpha*beta*xy[1])

  def phiY(self, xy, xyIm, g, alpha, beta, gamma2, gamma):
    return gamma * (1 + alpha*xy[0] + beta*xy[1]) * xyIm[1] - g * sqrt(gamma2 + beta**2) * xy[1]


class AtmosphereRT(Atmosphere):
  
  def __init__(self, g:float, alpha:float, beta:float, focal:float, error:float, nbPixels:int):
    super().__init(g, alpha, beta, focal, error, nbPixels)

        
