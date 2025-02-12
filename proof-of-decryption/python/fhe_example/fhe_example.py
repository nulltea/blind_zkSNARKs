import sys
sys.path.append('..')   # path to lazer module
from lazer import *     # import lazer python module
import hashlib          # for SHAKE128
import secrets          # for RNG
import numpy as np
from math import sqrt
from fhe_example_params import mod, deg, mod_t, v_e, ct_count

# BFV params
deg = 2048
# mod = 2**54+1 # CT modulus
# mod_t = 5 # Plaintext modulus (prime)

R = polyring_t(deg, mod)
Rt = polyring_t(deg, mod_t) # plaintext

# s = gen_binary_poly(size)
# to do 
s = poly_t.urandom_static(R, 2, secrets.token_bytes(32), 0)
print("BFV.key-gen.sk generated.")
# s.print()
a = poly_t.urandom_static(R, mod, secrets.token_bytes(32),0)
print("BFV.key-gen.a generated.")
# a.print()

# todo: std must be 3.19 (from Illia's repository)
e = poly_t.grandom_static(R,1, secrets.token_bytes(32), 0)
print("BFV.key-gen.e generated.")
# e.print()

b = -a*s + e
#b.print()
pk = (b, a)
print("BFV.key-gen.(pk,sk) generated.")

delta = round(mod/mod_t)
print("BFV.enc.delta generated.")
#print("delta = ",delta)

u = poly_t.urandom_static(R, 2, secrets.token_bytes(32), 0)
print("BFV.enc.u generated.")
# u.print()

e0 = poly_t.grandom_static(R,1, secrets.token_bytes(32), 0)
e1 = poly_t.grandom_static(R,1, secrets.token_bytes(32), 0)
print("BFV.enc.(e0,e1) generated.")

# -- encrypt
m = poly_t.urandom_static(Rt, 5, secrets.token_bytes(32), 0)
print("BFV.enc.m prepared randomly.")
# m.print()

# ct0 = poly_t(R, None)
# ct1 = poly_t(R, None)
# ct0.set_coeffs()
# ct0 = (pk[0]*u + e0 + delta*m) 
# ct1 = (pk[1]*u + e1). 

# print("m_delta=",m_delta.to_list()[:10])
m_modq = poly_t(R, m.to_list())
#print("m*delta=",m_modq.__mul__(delta).to_list()[:10])

m_delta = poly_t(R, (m_modq.__mul__(delta).to_list()))
#print("m*delta=",m_delta.to_list()[:10])
ct0 = (pk[0]*u + e0 + m_delta)
ct1 = (pk[1]*u + e1)

print("BFV.ct0 =",ct0.to_list()[:10])
print("BFV.ct1 =",ct1.to_list()[:10])
print("\n")

# --- decrypt
s_q = poly_t(R, (s.to_list()))
#print("s_q=",s_q.to_list()[:10])
ct_s = (ct0 + ct1*s_q)
ct_s_coeffs = ct_s.to_list()
#print("coeffs=",ct_s_coeffs[:10])

ct_s_coeffs = [coeff*mod_t for coeff in ct_s_coeffs]
#print("coeffs * t=",ct_s_coeffs[:10])
ct_s_coeffs = [round(coeff/(mod)) for coeff in ct_s_coeffs]
#print("coeffs / q=",ct_s_coeffs[:10])

#m_decrypted = [coeff % mod_t for coeff in ct_s_coeffs]
#print("coeffs mod t=",m_decrypted[:10])
m_decrypted = poly_t(Rt, ct_s_coeffs)

print("m=",m.to_list()[:10])
print("m_decrypted=",m_decrypted.to_list()[:10])

print(bool(m == m_decrypted))


# b = polyadd(polymul_ntt(negate(a), s, modulus, primitive_root, poly_mod), negate(e), modulus, poly_mod)
# return (b, a), s


# --------------- Proving R1 ----------------

# # B_e
# sigma_err = 1.55*2 # this set in grandom_static since we used log2o=1
# B_e = sqrt(2*deg)*sigma_err
# print("B_e = ", B_e)

# # B_v
# delta_phi_m = deg # to-do check!!?
# inf_norm_t = mod_t
# B_v = mod / (2*delta_phi_m*inf_norm_t) - 0.5
# print("B_v = ", B_v)

# v_inh
v_inh = ct0 + ct1 * s - m_delta

# e
e = pk[1]*s + pk[0]


#[-p1 , 1 , 0]
#[-c1 , 0 , 1]
# *
#[s , e , v]
# =
# p0 
# c0 - delta_m


A1 = polymat_t(R, 2, 1, None)
A1.set_col(0, polyvec_t(R, 2, [-pk[1], -ct1]))
A2 = polymat_t.identity(R, 2)
A = polymat_t(R, 2, 3, [A1, A2])

seed = b'\0' * 32
# w = polyvec_t(R, 3, [s,e,v_inh])
# w = polyvec_t(R, 3, [,s,s])
w = polyvec_t(R, 3)
w.brandom(1, seed, 0)

t = -A*w

# shake128 = hashlib.shake_128(bytes.fromhex("01"))
# P1PP = shake128.digest(32)      # proof system public randomness


from _fhe_example_params_cffi import lib 
prover = lin_prover_state_t(seed, lib.get_params("param"))
verifier = lin_verifier_state_t(seed, lib.get_params("param"))

# prover.set_statement(A, t)
# prover.set_witness(w)

# print("generate proof ...")
# proof = prover.prove()










# -------------------- vdec_ct.h generator --------------------

def add_header_guards(file_path="vdec_ct.h", guard_name="VDEC_CT_H"):
    header_guard = f"#ifndef {guard_name}\n#define {guard_name}\n#include <stdint.h>\n\n"
    try:
        with open(file_path, "x") as file:
            file.write(header_guard)
    except FileExistsError:
        pass

def add_footer_guards(file_path="vdec_ct.h", guard_name="VDEC_CT_H"):
    footer_guard = f"\n#endif /* {guard_name} */\n"
    with open(file_path, "a") as file:
        file.write(footer_guard)

def write_c_array_to_file(array, ct_count=1, file_path="vdec_ct.h", array_name="static_sk", array_type="uint8_t"):
    add_header_guards(file_path)
    
    # Convert the array elements to strings and join them with commas
    array_elements = ", ".join(map(str, array))
    array_length = len(array)
    
    # Read the file contents if it exists
    try:
        with open(file_path, "r") as file:
            existing_content = file.read()
    except FileNotFoundError:
        existing_content = ""

    # Check if the array definition already exists in the file
    if f"static const {array_type} {array_name}[]" in existing_content:
        # Append to the existing array definition
        existing_content = existing_content.rstrip().rstrip("};\n")
        array_elements = f", {array_elements}\n}};\n"
    else:
        # Create a new array definition
        if ct_count == 1:
            array_definition = f"static const {array_type} {array_name}[] = {{\n    {array_elements}\n}};\n"
        else:
            # Build an array of arrays
            array_of_arrays = ",\n    ".join([f"{{ {array_elements} }}" for _ in range(ct_count)])
            array_definition = f"static const {array_type} {array_name}[{ct_count}][{array_length}] = {{\n    {array_of_arrays}\n}};\n"
        array_elements = array_definition
    
    # Write or append the formatted array definition to the file
    with open(file_path, "a") as file:
        file.write(array_elements)



# write_c_array_to_file(s.to_list(), array_name="static_sk", array_type="int64_t")
# write_c_array_to_file(ct0.to_list(), array_name="static_ct0", array_type="int64_t")
# write_c_array_to_file(ct1.to_list(), array_name="static_ct1", array_type="int64_t")
# write_c_array_to_file(m_delta.to_list(), array_name="static_m_delta", array_type="int64_t")
# write_c_array_to_file(v_inh.to_list(), array_name="static_v_inh", array_type="int64_t")
# write_c_array_to_file(pk[0].to_list(), array_name="static_pk0", array_type="int64_t")
# write_c_array_to_file(pk[1].to_list(), array_name="static_pk1", array_type="int64_t")
# write_c_array_to_file(e.to_list(), array_name="static_e", array_type="int64_t")
# add_footer_guards()

write_c_array_to_file(s.to_list(), array_name="static_sk", array_type="int64_t")
write_c_array_to_file(ct0.to_list(), ct_count=ct_count, array_name="static_ct0", array_type="int64_t")
write_c_array_to_file(ct1.to_list(), ct_count=ct_count, array_name="static_ct1", array_type="int64_t")
write_c_array_to_file(m_delta.to_list(), ct_count=ct_count, array_name="static_m_delta", array_type="int64_t")
write_c_array_to_file(v_inh.to_list(), ct_count=ct_count, array_name="static_v_inh", array_type="int64_t")
write_c_array_to_file(pk[0].to_list(), array_name="static_pk0", array_type="int64_t")
write_c_array_to_file(pk[1].to_list(), array_name="static_pk1", array_type="int64_t")
write_c_array_to_file(e.to_list(), ct_count=ct_count, array_name="static_e", array_type="int64_t")
add_footer_guards()

