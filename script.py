from main import Atmosphere
import sys

# print("max rec", sys.getrecursionlimit())
sys.setrecursionlimit(2000)

atm = AtmosphereRT(g=10, alpha=0.01, beta=0.1, focal=16, error=1000, nbPixels=2000)
# print(atm)
atm.drawScreen("diapo")
atm.drawScreen("screen")
atm.drawScreen("target")
atm.drawScreen("corrected")
atm.drawScreen("final")
atm.prettyPrint()