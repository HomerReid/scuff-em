from numpy import *
from sys import argv

if __name__ == '__main__':
	xC, yC, zC = (0.1835, 0.4225, 0)

	dx = 0.01 # in um

	with open("EPFile", 'w') as fp:
		for x in arange(0, .1, dx):
			print("%.4f %.4f %.4f"%(x+xC, yC, zC), file=fp)