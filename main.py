from numpy import arctan, random
import numpy
from math import pi, cos, sin, sqrt
from matplotlib.pyplot import plot, show, title, xlabel, ylabel, subplot

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
    if not self.diapo: self.diapo = self.createDiapo()
    if not self.screen: self.screen = self.createScreen(self.diapo)
    # self.computeParameters()
    # if not self.screenTarget: self.screenTarget = self.createScreenTarget()
    # if not self.diapoCorrected: self.diapoCorrected = self.createDiapoCorrected()
    # if not self.screenFinal: self.screenFinal = self.createScreen(self.diapoCorrected, True)

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
    data = [int(value * 100 / nbPixels) for value in arrayPixels]
    names = []
    while lower < 1.11:
      names.append(f"Supérieures à {round(lower, 1)}")
      lower += 0.1
    ax = subplot(111)
    width=0.3
    bins = map(lambda x: x-width/2,range(1,len(data)+1))
    ax.bar(bins,data,width=width)
    ax.set_xticks(map(lambda x: x, range(1,len(data)+1)))
    ax.set_xticklabels(names,rotation=45, rotation_mode="anchor", ha="right")
    show()
    


  def computePrecision(self, value, approx):
    if value - approx != 0:
      if value != 0: return int(abs(value / (value - approx)))
      return int(abs(approx / (value - approx)))
    return "infinity"

  def drawScreen (self, plots, titleGraph="Pas de titre"):
      xs, ys = [], []
      for index in range(self.nbPoints):
        xs += [p[0] for p in plots[index]]
        ys += [p[1] for p in plots[index]]
      plot(xs, ys, 'o', color='black')
      title(titleGraph)
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
      diapo = []
    for indexY in range(self.nbPoints):
      newColumn, valueX = [], startX
      diapo.append(newColumn)
      for indexX in range(self.nbPoints):
        newColumn.append((valueX, valueY))
        valueX += stepX
      valueY += stepY
    return diapo

  def createScreen (self, diapo, error=True):
    cst = self.g/self.gamma
    numY = sqrt(self.beta**2 + self.gamma2)
    screenFinal = []
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
      screenFinal.append(screenLine)
    return screenFinal


  def createScreenTarget(self):
    startX = - self.gGd * self.sizeX / (2 * self.focal)
    stepX = -startX / ((self.nbPoints - 1) / 2)
    startY = - self.gGd * self.sizeY / (2 * self.focal)
    stepY = -startY / ((self.nbPoints - 1) / 2)
    valueY = startY
    screenTarget = []
    for indexY in range(self.nbPoints):
      newColumn, valueX = [], startX
      screenTarget.append(newColumn)
      for indexX in range(self.nbPoints):
        newColumn.append((valueX, valueY))
        valueX += stepX
      valueY += stepY
    return screenTarget

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
    diapoCorrected = []
    for line in self.screenTarget:
      screenLine = []
      for xy in line:
        numX = nxx * xy[0] + nxy * xy[1]
        den = dc + dx * xy[0] + dy * xy[1]
        numY = nyy * xy[1]
        screenLine.append((numX / den, numY / den))
      diapoCorrected.append(screenLine)
    return diapoCorrected

  def computeParameters(self):
    Atmosphere.gGd, Atmosphere.alphaGd, Atmosphere.betaGd = self.gradientDescent(self.costFunction, [0,0,0], 1)

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
    g,alpha,beta= array
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
  
  def __init__(self, g:float, alpha:float, beta:float, focal:float, error:float, nbPixels:int, rotationAxis:list=[0.0,0.0,1.0], rotationAngle:float=0.0):
    self.rotationAxis = rotationAxis
    self.rotationAngle = rotationAngle * pi / 180
    super().__init__(g, alpha, beta, focal, error, nbPixels)
    self.diapoRotate = self.createDiapoRotate(self.diapo, self.rotationAngle, self.rotationAxis)
    self.screenRotate = self.createScreen(self.diapoRotate)
    self.computeParameters()
    if not self.screenTarget: self.screenTarget = self.createScreenTarget()
    if not self.diapoCorrected: self.diapoCorrected = self.createDiapoCorrected()
    forScreenFinal = self.createDiapoRotate(self.diapoCorrected, self.rotationAngle, self.rotationAxis) if self.rotationAxis else self.diapoCorrected
    if not self.screenFinal: self.screenFinal = self.createScreen(forScreenFinal, False)


  def createDiapoRotate(self, diapo, angle, axis):
    diapoRotate = []
    for column in diapo:
      newColumn = []
      for (x,y) in column:
        newColumn.append(self.__rotation(x, y, angle, axis))
      diapoRotate.append(newColumn)
    return diapoRotate

  def createDiapoCorrected(self):
    alphaBetaCorrection = super().createDiapoCorrected()
    if not self.rotationAngle: return alphaBetaCorrection
    return self.createDiapoRotate(alphaBetaCorrection, -self.rotationAngleGd, self.rotationAxisGd)


  def __rotation(self, x, y, angle, axis):
    cosV = cos(angle)
    sinV = sin(angle)
    ux, uy, uz = axis
    xNew = (ux**2 * (1 - cosV) + cosV) * x + (ux*uy * (1 - cosV) - uz * sinV) * y + (ux*uz * (1 - cosV) + uy * sinV)
    yNew = (ux*uy * (1 - cosV)  + uz * sinV) * x + (uy**2 * (1 - cosV) + cosV) * y + (uy*uz * (1 - cosV) - ux * sinV)
    zNew = (ux*uz * (1 - cosV)  - uy * sinV) * x + (uy*uz * (1 - cosV) + ux * sinV) * y + (uz**2 * (1 - cosV) + cosV)
    return xNew / zNew, yNew / zNew

  def computeParameters(self):
    if not self.rotationAngle:
      Atmosphere.gGd, Atmosphere.alphaGd, Atmosphere.betaGd = self.gradientDescent(self.costFunction, [0,0,0], 1)
    else:
      params = self.gradientDescent(self.costFunction, [0,0,0,0,0,0], 1)
      AtmosphereRT.gGd, AtmosphereRT.alphaGd, AtmosphereRT.betaGd = params[0], params[1], params[2]
      AtmosphereRT.rotationAxisGd = [params[3], params[4], sqrt(1 - params[3]**2 - params[4]**2)]
      AtmosphereRT.rotationAngleGd = params[5]
      print("params", params, params[5] * 180 / pi)

  def costFunction (self, array):
    if not self.rotationAngle:
      return super().costFunction(array)
    g,alpha,beta, rotationAxisX, rotationAxisY, rotationAngle = array
    gamma2 = 1 + alpha**2
    gamma = sqrt(gamma2)
    cost = 0
    axis = [rotationAxisX, rotationAxisY, sqrt(1 - rotationAxisX**2 - rotationAxisY**2)]
    diapo = self.createDiapoRotate(self.diapo, rotationAngle, axis)
    for i in range(len(diapo)):
      line = diapo[i]
      for j in range(len(line)):
        phiX = self.phiX(line[j], self.screenRotate[i][j], g, alpha, beta, gamma2, gamma)
        phiY = self.phiY(line[j], self.screenRotate[i][j], g, alpha, beta, gamma2, gamma)
        cost += phiX**2 + phiY**2
    return cost

  def prettyPrint(self):
    super().prettyPrint()
    realAngle = self.rotationAngle * 180 / pi
    computedAngle = self.rotationAngleGd * 180 / pi
    print(f"precision angle : {self.computePrecision(realAngle, computedAngle)} precision axis X : {self.computePrecision(self.rotationAxis[0], self.rotationAxisGd[0])} precision axis Y : {self.computePrecision(self.rotationAxis[1], self.rotationAxisGd[1])}")



        
