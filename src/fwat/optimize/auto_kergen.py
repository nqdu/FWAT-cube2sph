import sympy as sp
import sys 

def wrap_string(s, max_len=60,indent=4):
    import re 
    tokens = re.split(r'(\*\*|[+\-*/=])', s)
    result = ''
    line = ''
    prefix = ''

    for i, token in enumerate(tokens):
        if len(line + token) > max_len:
            if line.strip():
                result += prefix + line.rstrip() + ' \\\n'
                prefix = ' ' * indent + ' '  # indent with '=' at new line
            line = token
        else:
            line += token

    result += prefix + line.strip()  # final line

    return result

def _c66mat(A,C,L,N,F,gc,gs):
    C0 = sp.Array([
        [A,A-2*N,F,0,0,0],
        [A-2*N,A,F,0,0,0],
        [F,F,C,0,0,0],
        [0,0,0,L - gc,-gs,0],
        [0,0,0,-gs,L+gc,0],
        [0,0,0,0,0,N]
    ])

    return C0

def dhti_call():
    vp,vs,rho,gsp,gcp = sp.symbols("vp,vs,rho,gsp,gcp")
    A = rho*vp**2 
    C = A 
    N = rho*vs**2
    L = N 
    gs = gsp * L 
    gc =  gcp * L 
    eta = 1
    F = eta * (A - 2 * L)
    C0 = _c66mat(A,C,L,N,F,gc,gs)

    # compute derivative 
    for ii,param in enumerate([vp,vs,rho,gcp,gsp]):
        expr = ""
        out = sp.diff(rho,param)
        if out !=0:
            expr += f"md_kl[-1,...]* ({str(sp.pycode(out))})+" 
        kk = 0
        for i in range(6):
            for j in range(i,6):
                out = sp.diff(C0[i,j],param)
                if out !=0:
                    expr = expr + f"md_kl[{kk},...]*" + f"({str(sp.pycode(out))})+"
                kk += 1
        if expr[-1] == "+":
            expr = expr[:len(expr)-1]
        print(f"md_newkl[{ii},...] = {wrap_string(expr)}")

def dtti_call():
    rho,vph,vpv,vsh,vsv,eta,gcp,gsp = sp.symbols("rho,vph,vpv,vsh,vsv,eta,gcp,gsp")
    A = rho * vph**2 
    L = rho * vsv**2 
    C = rho * vpv**2 
    N = rho * vsh**2
    gc = gcp * L 
    gs = gsp * L 
    F = eta * (A - 2* L)
    C0 = _c66mat(A,C,L,N,F,gc,gs)

    # compute derivative 
    for ii,param in enumerate([vph,vpv,vsh,vsv,rho,eta,gcp,gsp]):
        expr = ""
        out = sp.diff(rho,param)
        if out !=0:
            expr += f"md_kl[-1,...]* ({str(sp.pycode(out))})+" 
        kk = 0
        for i in range(6):
            for j in range(i,6):
                out = sp.diff(C0[i,j],param)
                if out !=0:
                    expr = expr + f"md_kl[{kk},...]*" + f"({str(sp.pycode(out))})+"
                kk += 1
        if expr[-1] == "+":
            expr = expr[:len(expr)-1]
        print(f"md_newkl[{ii},...] = {wrap_string(expr)}")

def dtti_call_with_kappa():
    rho,vp,vs,kappaa,kappab,eta,gcp,gsp = sp.symbols("rho,vp,vs,kappaa,kappab,eta,gcp,gsp")

    # compute vph,vpv,vsh,vsv
    # vp = 1/2 (vph + vpv)
    # kappaa = vph/vpv 
    # => vpv = vp * 2 /(1 + kappaa)
    # => vph = vpv * kappaa
    vpv = vp * 2 /(1 + kappaa)
    vsv = vs * 2 /(1 + kappab)
    vph = vpv * kappaa
    vsh = vsv * kappab
    A = rho * vph**2 
    L = rho * vsv**2 
    C = rho * vpv**2 
    N = rho * vsh**2
    gc = gcp * L 
    gs = gsp * L 
    F = eta * (A - 2* L)
    C0 = _c66mat(A,C,L,N,F,gc,gs)

    # compute derivative 
    for ii,param in enumerate([vp,vs,rho,kappaa,kappab,eta,gcp,gsp]):
        expr = ""
        out = sp.diff(rho,param)
        if out !=0:
            expr += f"md_kl[-1,...]* ({str(sp.pycode(out))})+" 
        kk = 0
        for i in range(6):
            for j in range(i,6):
                out = sp.diff(C0[i,j],param)
                if out !=0:
                    expr = expr + f"md_kl[{kk},...]*" + f"({str(sp.pycode(out))})+"
                kk += 1
        if expr[-1] == "+":
            expr = expr[:len(expr)-1]
        print(f"md_newkl[{ii},...] = {wrap_string(expr)}")

def dtti_call_with_voigt():
    rho,vp,vs,kappaa,kappab,eta,gcp,gsp = sp.symbols("rho,vp,vs,kappaa,kappab,eta,gcp,gsp")

    # compute vph,vpv,vsh,vsv
    # vp = sqrt ((vph**2*2 + vpv**2)/3)
    # kappaa = (vph - vpv)/vpv 
    # => vpv = 
    # => vph = vpv * kappaa
    vpv_sq = 3 * vp**2 / (2 * (1+kappaa)**2 + 1)
    vph_sq = (1 + kappaa)**2 * vpv_sq
    vsv_sq = 3 * vs**2 / (2 * (1+kappab)**2 + 1)
    vsh_sq = (1 + kappab)**2 * vsv_sq
    A = rho * vph_sq
    L = rho * vsv_sq 
    C = rho * vpv_sq 
    N = rho * vsh_sq 
    gc = gcp * L 
    gs = gsp * L 
    F = eta * (A - 2* L)
    C0 = _c66mat(A,C,L,N,F,gc,gs)

    # compute derivative 
    for ii,param in enumerate([vp,vs,rho,kappaa,kappab,eta,gcp,gsp]):
        expr = ""
        out = sp.diff(rho,param)
        if out !=0:
            expr += f"md_kl[-1,...]* ({str(sp.pycode(out))})+" 
        kk = 0
        for i in range(6):
            for j in range(i,6):
                out = sp.diff(C0[i,j],param)
                if out !=0:
                    expr = expr + f"md_kl[{kk},...]*" + f"({str(sp.pycode(out))})+"
                kk += 1
        if expr[-1] == "+":
            expr = expr[:len(expr)-1]
        print(f"md_newkl[{ii},...] = {wrap_string(expr)}")

def main():
    if len(sys.argv) != 3:
        print("Usage: ./this mdtype kltype")
        print("Example: ./this dtti 2")
        exit(1)
    
    mdtype = sys.argv[1]
    kltype = int(sys.argv[2])
    if (mdtype not in ['dtti']) or (kltype not in [1,2,3]):
        print(f"mdtype = {mdtype} and kltype = {kltype} is not implemented")

    if mdtype == "dtti":
        if kltype == 1:
            dhti_call()
        elif kltype == 2:
            dtti_call()
        elif kltype ==3:
            dtti_call_with_voigt()
        else:
            print(f"mdtype = {mdtype} and kltype = {kltype} is not implemented")


if __name__ == "__main__":
    main()