import sys
#bit_length = 64
#r = 3
def main():
	try:
		x_min = x_min = (bit_length / 2) - 20
		if (x_min < 0):
			x_min = 0
		x_max = bit_length / 2 + 20
		if (bit_length < x_max):
			x_max = bit_length
		for x in range (x_min, x_max):
			y_min = ceil((bit_length - x - 20)/log(r,2))
			if (y_min < 0):
				y_min = 0
			for y in range (y_min, floor((bit_length - x)/log(r,2))):
				for f in range (1, bit_length - x - ceil(y*log(r,2)), 2):
					if (f % r != 0 and is_prime(f * 2**x * r**y - 1)):
						print(f,x,y)
		print("")
	except:
		print("Please define variables:\nbit_length\nr")
main()