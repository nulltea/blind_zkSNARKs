import sys
sys.path.append('..')   # path to lazer module
from lazer import *     # import lazer python module
from lazer import _invmod
import hashlib      # for SHAKE128
import time
from labrador import *

FALCON_RING=polyring_t(512,12289)
BIGFALCON_RING=polyring_t(512,LAB_RING_32.mod)
BIGMOD_RING=polyring_t(256,LAB_RING_32.mod)
PRIMESIZE=str(math.ceil(math.log2(BIGMOD_RING.mod)))

# write the public keys and the v polynomials as a0+a1*FALC_SPLIT_BASE
FALC_SPLIT_BASE=2**7

# the falcon secret polynomials gets mapped to half the degree 
FALC_SEC_SPLIT_NORM=17017363//2

# the v polynomial gets mapped to half the degree and written in base FALC_SPLIT_BASE
V_SPLIT_NORM=round(1248245003*.75)//(2*FALC_SPLIT_BASE**2) 

PK_SPLIT_NORM=256*FALC_SPLIT_BASE**2

SIG_NORM_LIST=[PK_SPLIT_NORM]*4
SIG_NORM_LIST.extend([FALC_SEC_SPLIT_NORM]*4)
SIG_NORM_LIST.extend([PK_SPLIT_NORM]*2)
SIG_NORM_LIST.extend([V_SPLIT_NORM]*2)
print(SIG_NORM_LIST)

def center_mod(a,m):
    a=a % m
    if a>m//2:
        a=a-m
    return a

def decompose(pol:poly_t,base,loops=0):
    #pol.redc()
    pol.redp()
    if loops==0:
        loops=math.ceil(math.log(pol.ring.mod,base))
    temp_pol=poly_t(pol.ring)
    #cur=ffi.new("int64_t []",pol.ring.deg)
    #top=ffi.new("int64_t []",pol.ring.deg)
    #lib.poly_get_coeffvec_i64(top, pol.ptr)
    top=pol.make_i64array()
    res=polyvec_t(pol.ring,loops)
    for i in range(loops):
        top,cur=armod(top,pol.ring.deg,base)
        #lib.poly_set_coeffvec_i64(temp_pol.ptr,cur)
        temp_pol.set_i64array(cur)
        res.set_elem(temp_pol,i)
    return res

def recompose(polvec:polyvec_t,base):
    pol=poly_t(polvec.ring)
    curmul=1
    for i in range(polvec.dim):
        pol+=curmul*polvec.get_elem(i)
        curmul*=base
    return pol

def armod(vec_in,deg,mod,center=False):
    vec_out=ffi.new("int64_t []",deg)
    top=ffi.new("int64_t []",deg)
    for i in range(deg):
        # TODO switch to shifts later if mod is always a power of 2
        vec_out[i]=vec_in[i] % mod
        if center:
            vec_out[i]=center_mod(vec_in[i],mod)
        top[i]=(vec_in[i]-vec_out[i]) // mod
    return top,vec_out

def makeGvec(ring,base,dim):
    G=polyvec_t(ring,dim)
    for i in range(dim):
        G.set_elem(poly_t(ring,{0:base**i}),i)
    return G

def make_LRhash(ring:polyring_t,length,seed):
    Lhash=polyvec_t.urandom_bnd_static(ring,length,0,ring.mod-1,seed,0)
    Rhash=polyvec_t.urandom_bnd_static(ring,length,0,ring.mod-1,seed,1)
    return Lhash,Rhash

def check_path(root:poly_t, node:poly_t, i:int, path:list, dec_base: int, Lhash:polyvec_t, Rhash: polyvec_t):
    temp=poly_t(root.ring)
    count=0
    hash_length=Lhash.dim
    while i != 0:
        if i%2 == 0:
            node=Lhash*decompose(path[count],dec_base,hash_length) + \
                 Rhash*decompose(node,dec_base,hash_length)
        else:
            node=Lhash*decompose(node,dec_base,hash_length) + \
                 Rhash*decompose(path[count],dec_base,hash_length)
        count+=1
        i=(i-1)//2
    return node==root

def make_falcon_pk_leaves(num_leaves:int,Lhash:polyvec_t):
    falcon_pk=[]
    falcon_sig=[]
    leaves=[]
    base=FALC_SPLIT_BASE
    X=poly_t(BIGMOD_RING,{1:1})
    inv_fal_mod=_invmod(12289,BIGFALCON_RING.mod)
    
    shake128 = hashlib.shake_128(bytes.fromhex("44"))
    TARGPP = shake128.digest(32)
    f_t=poly_t.urandom_static(FALCON_RING,FALCON_RING.mod,TARGPP,0)
    l_t=f_t.lift(BIGFALCON_RING)


    i=0

    while i <num_leaves:
        skenc,pkenc,pkpol=falcon_keygen()
        l_s1,l_s2=falcon_preimage_sample(skenc,l_t)
        l_s1=l_s1.lift(BIGFALCON_RING)
        l_s2=l_s2.lift(BIGFALCON_RING)
        l_pk=pkpol.lift(BIGFALCON_RING)
        
        v=poly_t(BIGFALCON_RING)
        v=(l_t-l_s1-l_pk*l_s2)*inv_fal_mod
        var=v.make_i64array()
        var_high,var_low=armod(var,v.ring.deg,base)
        v_high=poly_t(BIGFALCON_RING)
        v_low=poly_t(BIGFALCON_RING)
        v_high.set_i64array(var_high)
        v_low.set_i64array(var_low)
        assert v_low+v_high*base == v
        assert v_low.linf()<base and v_high.linf()<base


        l_pkar=l_pk.make_i64array()
        pkar_high,pkar_low = armod(l_pkar,l_pk.ring.deg,base)
        l_pk_low=poly_t(BIGFALCON_RING)
        l_pk_high=poly_t(BIGFALCON_RING)
        l_pk_low.set_i64array(pkar_low)
        l_pk_high.set_i64array(pkar_high)
        assert l_pk_low+l_pk_high*base == l_pk
        assert (l_pk_low+l_pk_high*base)*l_s2+l_s1 + 12289*(v_low+v_high*base) == l_t 
        assert l_pk_low.linf() < base and l_pk_high.linf()<base

        pk_split_low=l_pk_low.to_isoring(BIGMOD_RING)
        pk_split_high=l_pk_high.to_isoring(BIGMOD_RING)
        v_split_low=v_low.to_isoring(BIGMOD_RING)
        v_split_high=v_high.to_isoring(BIGMOD_RING)
        s1_split=l_s1.to_isoring(BIGMOD_RING)
        s2_split=l_s2.to_isoring(BIGMOD_RING)
        t_split=l_t.to_isoring(BIGMOD_RING)
        
        assert pk_split_low.l2sqr() < PK_SPLIT_NORM and pk_split_high.l2sqr() and v_split_low.l2sqr() < PK_SPLIT_NORM
        if v_split_high.get_elem(0).l2sq() > V_SPLIT_NORM or v_split_high.get_elem(1).l2sq() > V_SPLIT_NORM :
            print("V TOO BIG")
            continue
        if s1_split.get_elem(0).l2sq() > FALC_SEC_SPLIT_NORM or s1_split.get_elem(1).l2sq() > FALC_SEC_SPLIT_NORM \
            or s2_split.get_elem(0).l2sq() > FALC_SEC_SPLIT_NORM or s2_split.get_elem(1).l2sq() > FALC_SEC_SPLIT_NORM:
            print("S TOO BIG")
            continue
       
        assert (pk_split_low.get_elem(0)+base*pk_split_high.get_elem(0))*s2_split.get_elem(0) + \
            (pk_split_low.get_elem(1)+base*pk_split_high.get_elem(1))*s2_split.get_elem(1)*X + \
            s1_split.get_elem(0)+12289*(v_split_low.get_elem(0)+v_split_high.get_elem(0)*base) == t_split.get_elem(0)

        assert (pk_split_low.get_elem(1)+base*pk_split_high.get_elem(1))*s2_split.get_elem(0) + \
            (pk_split_low.get_elem(0)+base*pk_split_high.get_elem(0))*s2_split.get_elem(1) + \
            s1_split.get_elem(1)+12289*(v_split_low.get_elem(1)+v_split_high.get_elem(1)*base) == t_split.get_elem(1)

        falcon_pk.append(polyvec_t(BIGMOD_RING,4,[pk_split_low,pk_split_high]))
        falcon_sig.append(polyvec_t(BIGMOD_RING,8,[s1_split,s2_split,v_split_low,v_split_high]))
        assert Lhash.dim==4
        leaves.append(Lhash*falcon_pk[i])
        i+=1
    return leaves,falcon_pk,falcon_sig



class hash_tree:
    def __init__(self,ring:polyring_t,depth:int,leaves:list,dec_base:int,seed,LRhash=None):
        zpol=poly_t(ring)
        self.tree=[zpol]*(2**(depth+1)-1) # make an empty tree
        self.ring=ring
        self.depth=depth
        self.dec_base=dec_base
        self.hash_length=math.ceil(math.log(ring.mod,dec_base))
        if LRhash == None:
            self.Lhash,self.Rhash = make_LRhash(ring,self.hash_length,seed)
        else:
            self.Lhash,self.Rhash = LRhash
        self.leaves=leaves
        self.leaf_start=2**depth-1
        for i in range(len(leaves)):
            self.tree[self.leaf_start+i]=leaves[i]
        for i in range(2**depth-2,-1,-1):
            self.tree[i]=self.Lhash*decompose(self.tree[2*i+1],dec_base,self.hash_length) + \
                         self.Rhash*decompose(self.tree[2*i+2],dec_base,self.hash_length)
        
    def get_path(self,i):
        path=[]
        while i != 0:
            if i%2==0:
                path.append(self.tree[i-1])
            else:
                path.append(self.tree[i+1])
            i=(i-1)//2
        return path
    
    def decomposed_compath(self,i):
        path=self.get_path(i)
        node=poly_t(self.ring,self.tree[i])
        compath=[]
        pospath=[]
        count=0
        while i != 0:
            if i%2 == 0:
                Ldec=decompose(path[count],self.dec_base,self.hash_length)
                Rdec=decompose(node,self.dec_base,self.hash_length)
                compath.append(Ldec)
                compath.append(Rdec)
                if i>2:
                    node=self.Lhash*Ldec + self.Rhash*Rdec
                if(count>0):
                    pospath.append(1)
            else:
                Ldec=decompose(node,self.dec_base,self.hash_length)
                Rdec=decompose(path[count],self.dec_base,self.hash_length)
                compath.append(Ldec)
                compath.append(Rdec)
                if i>2:
                    node=self.Lhash*Ldec + self.Rhash*Rdec
                if(count>0):
                    pospath.append(0)
            i=(i-1)//2
            count+=1
        #for i in range(len(compath)):
        #    compath[i].redc()
        return compath,pospath

def test_proof():
    small_deg=64
    deg_list=[small_deg]
    num_pols_list=[2]
    norm_list=[2**19]
    num_constraints=1
    PS=proof_statement(deg_list,num_pols_list,norm_list,num_constraints,PRIMESIZE)
    R256=polyring_t(small_deg,LAB_RING_32.mod)
    R512=polyring_t(small_deg*2,LAB_RING_32.mod)
    shake128 = hashlib.shake_128(bytes.fromhex("00"))
    seed=shake128.digest(32)
    hash=polyvec_t.urandom_bnd_static(R256,2,0,10,seed,0)
    pol=poly_t.urandom_static(R512,10,seed,0)
    pol256=polyvec_t(R256,2)
    pol256=pol.to_isoring(R256) # does the same thing as the (commented out) code below
    for i in range(R256.deg):
        pol256.set_elem(pol.get_coeff(2*i),0,i)
        pol256.set_elem(pol.get_coeff(2*i+1),1,i)
    PS.fresh_statement([hash],[pol],hash*pol256)
    print(pol)
    print(pol256.get_elem(0))
    print(pol256.get_elem(1))
    PS.smpl_verify()


def create_proof(HT:hash_tree,node_list:list,falcon_pk=[],falcon_sig=[]):
    FALC=len(falcon_pk) > 0 and len (falcon_sig) > 0
    LN=len(node_list)
    for i in range(LN):
        assert node_list[i] > 2**HT.depth-2 and node_list[i]<2**(HT.depth+1)-1
    deg_list=[HT.ring.deg]*(12*FALC+2*HT.depth)*LN
    num_pols_list=[HT.hash_length]*2*HT.depth*LN
    
    norm_list=[math.ceil(HT.ring.deg*HT.hash_length*(HT.dec_base**2))]*2*HT.depth*LN
    num_constraints=HT.depth*LN
    if FALC:
        num_pols_list.extend([1]*12*LN) # 12 polynomials per falcon signature
        norm_list.extend(SIG_NORM_LIST*LN) # 12 signature polynmials
        #num_constraints+=(3*LN) # 1 for pk hashing and 2 for signature
        num_constraints+=LN # 1 for pk hashing. no signature yet 

    print(len(deg_list)," ",len(norm_list)," ",len(num_pols_list))
    PS=proof_statement(deg_list,num_pols_list,norm_list,num_constraints,PRIMESIZE)
    polzero=poly_t(HT.ring)

    negG=-makeGvec(HT.ring,HT.dec_base,HT.hash_length)
    negG.redc()
    cur_start=0
    for i in range(LN):
        print(PS.cur_witness_num)
        compath,pospath=HT.decomposed_compath(node_list[i])
        for j in range(len(compath)):
            PS.append_witness(compath[j])
        count=0
        print("length of compath=",len(compath))
        for j in range(0,len(compath),2):
            if j+2 < len(compath):
                stat_left=[HT.Lhash,HT.Rhash,negG]
                stat_right=polzero
                wit=[cur_start+j,cur_start+j+1,cur_start+j+2+pospath[count]]
            else:
                stat_left=[HT.Lhash,HT.Rhash]
                stat_right=HT.tree[0]
                wit=[cur_start+j,cur_start+j+1]
            
            PS.append_statement(stat_left,wit,stat_right)
            count+=1
        cur_start+=len(compath)
    
    START=2*HT.depth*LN
    if FALC:
        for i in range(LN):
            for j in range(4):
                PS.append_witness(falcon_pk[i].get_elem(j))
            stat_left=[negG,HT.Lhash.get_elem(0),HT.Lhash.get_elem(1),HT.Lhash.get_elem(2),HT.Lhash.get_elem(3)]
            #stat_right=HT.tree[node_list[i]]
            stat_right=polzero
            offset = (node_list[i] % 2 == 0)
            wit=[2*HT.depth*i+offset,START+12*i,START+12*i+1,START+12*i+2,START+12*i+3]
            PS.append_statement(stat_left,wit,stat_right)
            for j in range(8):
                PS.append_witness(falcon_sig[i].get_elem(j))



    PS.smpl_verify()
    stmnt=PS.output_statement()
    proof = PS.pack_prove()
    pack_verify(proof,stmnt,PRIMESIZE)

def main():
    

    shake128 = hashlib.shake_128(bytes.fromhex("05"))
    seed=shake128.digest(32)

    lstart=time.perf_counter()
    base=2**8
    loops=math.ceil(math.log(BIGMOD_RING.mod,base))
    leaves=[]
    depth=10
    LRhash=make_LRhash(BIGMOD_RING,4,seed)
    leaves,falcon_pk,falcon_sig = make_falcon_pk_leaves(2**depth,LRhash[0])
    #for i in range(2**depth):
    #    leaves.append(poly_t.urandom_static(BIGMOD_RING,BIGMOD_RING.mod,seed,i))
    HT=hash_tree(BIGMOD_RING,depth,leaves,base,seed,LRhash)
    # for i in range(len(HT.tree)):
    #     print(i)
    #     print(HT.tree[i])
    ind=2**depth-1+3
    path=HT.get_path(ind)
    cp=check_path(HT.tree[0],HT.tree[ind],ind,path,base,HT.Lhash,HT.Rhash)
    print(cp)
    path_start=time.perf_counter()
    compath,pospath=HT.decomposed_compath(ind)
    path_end=time.perf_counter()
    # for i in range(len(compath)):
    #     print(i)
    #     compath[i].print()
    G=makeGvec(BIGMOD_RING,base,loops)
    count=0
    for i in range(0,len(compath),2):
        if i+2<len(compath):
            print("-----")
            #print(HT.Lhash*compath[i]+HT.Rhash*compath[i+1]-G*compath[i+2+pospath[count]])
            print(HT.Lhash*compath[i]+HT.Rhash*compath[i+1]==G*compath[i+2+pospath[count]])
            count+=1
        else:
            print("-----")
            print(HT.Lhash*compath[i]+HT.Rhash*compath[i+1]==HT.tree[0])
    
    proof_start=time.perf_counter()
    index_list=[ind,ind+1,ind+4,ind+5,ind+9,ind+10,ind,ind+1,ind+4,ind+5,ind+9,ind+10]
    pk_list=[]
    sig_list=[]
    for i in index_list:
        pk_list.append(falcon_pk[i-2**depth+1])
        sig_list.append(falcon_sig[i-2**depth+1])
    create_proof(HT,index_list,pk_list,sig_list)
    proof_end=time.perf_counter()
    # for i in range(len(path)):
    #     print(path[i])
    # for i in range(1000):
    #     pol=poly_t.urandom_static(LAB_RING,LAB_RING.mod,seed,i)
    #     res=decompose(pol,2**8-1)
    #     #pol2=recompose(res,2**8-1)
    #     #print(pol-pol2)
    # lend=time.perf_counter()

    print(path_end-path_start)
    print(proof_end-proof_start)
    #test_proof()
    #make_falcon_pk_leaves(20,HT.Lhash)

if __name__ == "__main__":
    main()