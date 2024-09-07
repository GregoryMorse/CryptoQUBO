import itertools, operator, functools

qubo_path = "QuboTests/"

def square(x): return x*x
def inrange(L, H, x): return L <= x and x <= H
def validate_xor():
  for n in range(1, 18):
    tmax = (n+1)//2
    nsubs = (tmax-1).bit_length() #floor(log2(tmax))
    for X in itertools.product(*(((0, 1),)*n)):
      sum_X = sum(X)
      assert functools.reduce(operator.xor, X) == any(
        square(sum_X-T) == 0 #always positive by integer squaring
        for S in itertools.product(*(((0, 1),)*nsubs))
        for T in (1 if nsubs==0 else #for n==1 and n==2
                  2*(1 + sum(S[i-1]<<(i-1) for i in range(1, nsubs+1)))-1 if tmax&(tmax-1) == 0 else #tmax is an exact power of 2
                  2*(1 + (tmax - (1<<(nsubs-1)))*S[nsubs-1] + sum(S[i-1]<<(i-1) for i in range(1, nsubs)))-1,) #general case
      )
def validate_range():
  for n in range(1, 12):
    for L in range(0, n+1):
      for H in range(L, n+1):
        tmax = (H-L+1+1)//2
        nsubs = (tmax-1).bit_length() #floor(log2(tmax))
        if L==0 and H==n: continue #always True case
        for X in itertools.product(*(((0, 1),)*n)):
          sum_X = sum(X)
          assert inrange(L, H, sum_X) == any(
            (sum_X-f1)*(sum_X-f2) == 0 #always positive by root squeezing theorem since f2-f1==1
            for S in itertools.product(*(((0, 1),)*nsubs))
            for f1 in (L if nsubs==0 else
                       L+(H-L+1-2*(1<<(nsubs-1)))*S[nsubs-1] + sum(S[i-1]<<i for i in range(1, nsubs)),)
            for f2 in (f1 if nsubs== 0 and (H-L) & 1==0 else f1+1,)
          )
#validate_xor()
#validate_range()
"""
from sat_qubo import *
from crypto import *
def gen_proof_hash(smallcoeffs=False):
  for hashinfo in (MD5_hashinfo, SHA1_hashinfo, SHA256_hashinfo):
    hashlen, hashfunc, hashname, oraclefunc = hashinfo
    message = var_to_bits('m', 64*8) #var_to_bits('m', 55*8)
    M = bytes("e=mc^2", "ascii")
    h = oraclefunc(bytearray(M))
    Mblock = M + b'\x80' + b'\0'*(55-len(M)) + (len(M)*8).to_bytes(8, 'little' if hashname == "md5" else 'big')
    secret = [z for x in Mblock for z in integer_to_bits(x, 8)]
    hbits = [z for x in bytes.fromhex(h) for z in integer_to_bits(x, 8)]
    cnf, varmap, submap, formmaps, bqm3, const3, varmapq, submapq, result = get_hash_instance_qubo(hashinfo, message, merkledamgard=False, hashbits=hbits, smallcoeffs=smallcoeffs)
    solmap = {x[0]: y for x, y in zip(message, secret) if not x in (True, False)}
    #solmap.update({x[0]: y for x, y in zip(result, hbits)})
    solmap = wff_submap_to_sol(solmap, submap)
    assert check_sat(cnf, {varmap[x] if solmap[x] else -varmap[x] for x in solmap if x in varmap})
    print(len(bqm3), const3, len(submapq), len(varmapq))
    solmapq = qubo_submap_to_sol(solmap, submapq)
    from groqubo.solver_utils.qubo_ising import get_Q_from_BQM, qubo_sol_tobitstring, qubo_energy_func_quadratic
    import numpy as np
    Q, h_qubo, variable_map = get_Q_from_BQM( bqm3 )
    Q = Q + np.diag(h_qubo)
    revmap = {variable_map[k]: k for k in variable_map}
    solution = np.array([solmapq[revmap[i]] for i in range(len(Q))], dtype=Q.dtype)
    print(Q.shape[0], np.sum(Q != 0), np.sum(Q != 0)/(Q.shape[0]*Q.shape[1]), np.max(np.abs(Q)))
    print(qubo_energy_func_quadratic(solution, Q, const3), np.sum(solution).astype(np.uint64))
    bqm_to_qubo_file(bqm3, variable_map, qubo_path + hashname + ("-smallcoeffs" if smallcoeffs else "") + "-pw-" + M.hex() + ".qubo",
                    [f"{hashname}(Password)=={h} Password {len(M)} characters",
                      "Password Solution: " + str(M, 'ascii'), "MIN problem, Optimal Energy: " + str(-const3),
                      f"Message Indices: {', '.join(str(variable_map[x[0]]) for x in message)}",
                      f"Password Values: bits_to_integer({[x if x in (False, True) else 'x['+str(variable_map[x[0]])+']' for x in message]}).to_bytes({len(M)}, 'little')",
                      "def bits_to_integer(b): return sum(int(x)*(1<<i) for i, x in enumerate(b))"])
    with open(qubo_path + hashname + ("-smallcoeffs" if smallcoeffs else "") + "-pw-" + M.hex() + ".sol", 'w') as f:
      f.write(qubo_sol_tobitstring(solution))
def gen_proof_aes():
    for keysize, is_enc in ((128, True), (128, False), (192, True), (192, False), (256, True), (256, False)):
        data, key = b'e=mc^2!!!!!!!!!\x00', b'abcdefghijklmnopabcdefghijklmnop'[:keysize//8]
        oracle = aes_encrypt(data, key) if is_enc else aes_decrypt(data, key)
        oraclebits = [y for x in oracle for y in integer_to_bits(x, 8)]
        bqm, const, cnf, varmap, submap, submapq, keysub, datasub, res = aes_qubo(keysize, is_enc, rescmp=oraclebits)        
        solmap = {k[0]: v for k, v in zip(keysub, [y for x in key for y in integer_to_bits(x, 8)])}
        solmap.update({k[0]: v for k, v in zip(datasub, [y for x in data for y in integer_to_bits(x, 8)])})
        #solmap.update({k[0]: v for k, v in zip(res, oraclebits)})
        solmap = wff_submap_to_sol(solmap, submap)
        assert check_sat(cnf, {varmap[x] if solmap[x] else -varmap[x] for x in solmap if x in varmap})
        solmap = qubo_submap_to_sol(solmap, submapq)
        from groqubo.solver_utils.qubo_ising import get_Q_from_BQM, qubo_sol_tobitstring, qubo_energy_func_quadratic
        import numpy as np
        Q, h_qubo, variable_map = get_Q_from_BQM( bqm )
        Q = Q + np.diag(h_qubo)
        revmap = {variable_map[k]: k for k in variable_map}
        solution = np.array([solmap[revmap[i]] for i in range(len(Q))], dtype=Q.dtype)
        print(Q.shape[0], np.sum(Q != 0), np.sum(Q != 0)/(Q.shape[0]*Q.shape[1]), np.max(np.abs(Q)))
        print(qubo_energy_func_quadratic(solution, Q, const), np.sum(solution).astype(np.uint64))
        bqm_to_qubo_file(bqm, variable_map, qubo_path + "AES" + str(keysize) + "-" + ("enc" if is_enc else "dec") + "-" + data.hex() + "-" + key.hex() + ".qubo",
                        [f"AES{keysize}({data}, key)=={oracle.hex()}",
                          "Key Solution: " + str(key, 'ascii'), "MIN problem, Optimal Energy: " + str(-const),
                          f"Message Indices: {', '.join(str(variable_map[x[0]]) for x in datasub)}",
                          f"Key Indices: {', '.join(str(variable_map[x[0]]) for x in keysub)}",
                          f"Key Values: bits_to_integer({[x if x in (False, True) else 'x['+str(variable_map[x[0]])+']' for x in keysub]}).to_bytes({keysize//8}, 'little')",
                          "def bits_to_integer(b): return sum(int(x)*(1<<i) for i, x in enumerate(b))"])
        with open(qubo_path + "AES" + str(keysize) + "-" + ("enc" if is_enc else "dec") + "-" + data.hex() + "-" + key.hex() + ".sol", 'w') as f:
          f.write(qubo_sol_tobitstring(solution))
#gen_proof_hash()
#gen_proof_hash(smallcoeffs=True)
#gen_proof_aes()
"""

import numpy as np
import re
def load_qubo_problem(filename, sign=1):
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
    return np.array([int(x) for x in s], dtype=np.float_)
def qubo_energy_func_quadratic(q, Q, C_qubo=0):
    return np.sum(q * (Q @ q), axis=0) + C_qubo
def integer_to_bits(n, l=None):
    return [True if (n & (1<<i)) != 0 else False for i in range(n.bit_length() if l is None else l)]
def validate_qubos():
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
validate_qubos()
