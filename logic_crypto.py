import itertools, operator, functools

# The path of the QUBO test files and solutions
qubo_path = "QuboTests/"


import numpy as np
import re



def load_qubo_problem(filename, sign=1):
	"""
	Call to load a QUBO matrix from a file
		 
       
	Parameters
	----------
	filename : string
		The path to the QUBO file

	sign : int (optional)
		Set sign to 1 (default) if the QUBO problem is a minimization problem, or to -1 for maximization problem (which is converted to minimization by the -1 factor)

	Return
	------
	Returns with the following packed values: 
		qubo : standing for the QUBO matrix,
		C : standing for a constant offset,
		indices : the indices of the QUBO variables paired with keyindices
		keyindices : the labels of the QUBO variables loaded from the file
		sz : the number of QUBO variables
		min_val : the minimal value in the QUBO matrix
		max_val : the maximal value in the QUBO matrix
	"""


	qubo = None
	C = 0
	sz = None
	max_val = 0
	min_val = 0
	keyindices = None

	with open(filename, 'r') as f:
		for line in f:
			if line.startswith('p'):
				sz = int(line.split()[3])
				qubo = np.zeros((sz,sz), dtype=np.float64)
			elif line.startswith('c'):
				m = re.match("(.*) problem, Optimal Energy: (.*)", line)
				if not m is None: C = -int(m.group(2))

				m = re.match(".* Message Indices: (.*)", line)
				if not m is None: indices = [int(x) for x in m.group(1).split(", ")]

				m = re.match(".* Key Indices: (.*)", line)
				if not m is None: keyindices = [int(x) for x in m.group(1).split(", ")]
				continue
			else:
				i, j, coeff = line.split()[:3]
				i, j, coeff = int(i), int(j), np.float64(coeff)

				if coeff > max_val:
					max_val = coeff
				elif coeff < min_val:
					min_val = coeff

				if i == j:
					qubo[i][j] = coeff * sign  #change sign depending on the input
				else:
					# qubo[i][j] = coeff * sign #change sign
					qubo[i][j] = coeff/2.0 * sign #change sign
					qubo[j][i] = coeff/2.0 * sign #change sign


	return qubo, C, indices, keyindices, sz, min_val, max_val




def qubo_bitstring_to_sol(s):
	"""
	Call to convert a bitstring input into np array of the solutions
		 
       
	Parameters
	----------
	s : string
		The bitstring of a potential QUBO vector


	Return
	------
	Returns with the QUBO variable values stored in an nunpy array.
	"""

	return np.array([int(x) for x in s], dtype=np.float_)



def qubo_energy_func_quadratic(q, Q, C_qubo=0):
	"""
	Call to evaluate the energy of the QUBO problem for the given QUBO vector q
		 
       
	Parameters
	----------
	q : numpy array
		The array containing the QUBO variables for which the energy is calculated

	Q : numpy array
		The QUBO matrix

	C_qubo : constant number (optional)
		The QUBO energy is offsetted by t his value


	Return
	------
	Returns with the calculated QUBO energy.
	"""

	return np.sum(q * (Q @ q), axis=0) + C_qubo



def integer_to_bits(n, l=None):
	"""
	Call to determine the bits of an integer and put them into a list
		 
       
	Parameters
	----------
	n : int
		The integer inout for which the bits are determined

	l : int (optional)
		The count of the tracked bits. if None, the bit-length of n is taken.


	Return
	------
	Returns with the list of bits.
	"""


	return [True if (n & (1<<i)) != 0 else False for i in range(n.bit_length() if l is None else l)]




def validate_qubos():
	"""
	Call to validate the example QUBO instances
		 
       

	Return
	------
	Returns with the list of bits.
	"""


	for qubofile in ("md5-pw-653d6d635e32", "sha1-pw-653d6d635e32", "sha256-pw-653d6d635e32",
                   "md5-smallcoeffs-pw-653d6d635e32", "sha1-smallcoeffs-pw-653d6d635e32", "sha256-smallcoeffs-pw-653d6d635e32",
                   "AES128-enc-653d6d635e3221212121212121212100-6162636465666768696a6b6c6d6e6f70",
                   "AES128-dec-653d6d635e3221212121212121212100-6162636465666768696a6b6c6d6e6f70",
                   "AES192-enc-653d6d635e3221212121212121212100-6162636465666768696a6b6c6d6e6f706162636465666768",
                   "AES192-dec-653d6d635e3221212121212121212100-6162636465666768696a6b6c6d6e6f706162636465666768",
                   "AES256-enc-653d6d635e3221212121212121212100-6162636465666768696a6b6c6d6e6f706162636465666768696a6b6c6d6e6f70",
                   "AES256-dec-653d6d635e3221212121212121212100-6162636465666768696a6b6c6d6e6f706162636465666768696a6b6c6d6e6f70"):

		Q, C, indices, keyindices, _, _, _ = load_qubo_problem(qubo_path + qubofile + ".qubo")

		with open(qubo_path + qubofile + ".sol") as f:
			sol = qubo_bitstring_to_sol(f.readline())

		if qubofile.startswith("AES"):
			keysize = int(qubofile[3:6])
			Mblock, key = b'e=mc^2!!!!!!!!!\x00', b'abcdefghijklmnopabcdefghijklmnop'[:keysize//8]
			keybits = [y for x in key for y in integer_to_bits(x, 8)]

			for idx, bit in zip(keyindices, keybits):
				assert int(sol[idx]) == bit

		else:
			M = bytes("e=mc^2", "ascii")
			Mblock = M + b'\x80' + b'\0'*(55-len(M)) + (len(M)*8).to_bytes(8, 'little' if qubofile.startswith("md5") else 'big')


		secret = [z for x in Mblock for z in integer_to_bits(x, 8)]

		for idx, bit in zip(indices, secret):
			assert int(sol[idx]) == bit

		print(qubofile, "Energy with correct message: ", qubo_energy_func_quadratic(sol, Q, C))


# validate the QUBO examples
validate_qubos()
