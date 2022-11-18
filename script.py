from main import AtmosphereRT
from math import sqrt
import sys

# print("max rec", sys.getrecursionlimit())
sys.setrecursionlimit(10000)

atm = AtmosphereRT(g=10, alpha=0.01, beta=0.1, focal=16, error=1000, nbPixels=2000, rotationAxis=[0.1,0.1,sqrt(0.98)], rotationAngle=10)
# atm = AtmosphereRT(g=10, alpha=0.01, beta=0.1, focal=16, error=1000, nbPixels=2000)
atm.drawScreen(atm.diapo, "diapositive initiale")
if atm.rotationAngle:
  atm.drawScreen(atm.diapoRotate, "diapositive après rotation")
  atm.drawScreen(atm.screenRotate, "écran après rotation")
atm.drawScreen(atm.screenTarget, "écran cible")
atm.drawScreen(atm.diapoCorrected, "diapositive corrigée")
atm.drawScreen(atm.screenFinal, "écran corrigée")
atm.prettyPrint()