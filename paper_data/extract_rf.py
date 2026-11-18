"""
extract_rf.py - provenance of the CSVs for Figs. 7(c), 8, 9, 10 and 11 (Steps 10b, 12a, 12b): reads the plotted positions
straight from the vector drawing of the paper's PDF [Nazir et al., IEEE TED 72(9) 2025] with PyMuPDF.
Run:  python3 paper_data/extract_rf.py PATH/TO/paper.pdf [OUTDIR]      (needs: pip install pymupdf)
Axis mapping: chart frames and tick marks found in the drawing (frame(), ticks_len(), right_ticks()); Smith and
polar charts by their outline circles. Markers are small closed Bezier shapes, model lines long polylines.
The PDF itself is not in the repository (copyright); figures 6(a)-7(b) were extracted earlier with the same helpers.
"""
import pymupdf as fitz, numpy as np, collections
import sys
D=fitz.open(sys.argv[1] if len(sys.argv) > 1 else 'paper.pdf')
def frame(pg, r):
    p=D[pg-1]; H=[];V=[]
    for dr in p.get_drawings():
        if not r.contains(dr["rect"]) or not dr["color"] or max(dr["color"])>0.05: continue
        for it in dr["items"]:
            if it[0]!='l': continue
            a,b=it[1],it[2]
            if abs(a.y-b.y)<0.01: H.append((abs(a.x-b.x),min(a.x,b.x),max(a.x,b.x),a.y))
            if abs(a.x-b.x)<0.01: V.append((abs(a.y-b.y),min(a.y,b.y),max(a.y,b.y),a.x))
    H.sort(reverse=True); V.sort(reverse=True)
    return H, V
def ticks(pg, r, fr, L=3.0):
    """short black tick lines touching the frame"""
    x0,x1,y0,y1=fr; p=D[pg-1]; xt=set(); yl=set(); yr=set()
    for dr in p.get_drawings():
        if not r.contains(dr["rect"]) or not dr["color"] or max(dr["color"])>0.05: continue
        for it in dr["items"]:
            if it[0]!='l': continue
            a,b=it[1],it[2]
            if abs(a.x-b.x)<0.01 and abs(a.y-b.y)<L and (abs(max(a.y,b.y)-y1)<0.05 or abs(min(a.y,b.y)-y1)<0.05): xt.add(round(a.x,2))
            if abs(a.y-b.y)<0.01 and abs(a.x-b.x)<L:
                if abs(min(a.x,b.x)-x0)<0.05: yl.add(round(a.y,2))
                if abs(max(a.x,b.x)-x1)<0.05: yr.add(round(a.y,2))
    return sorted(xt), sorted(yl), sorted(yr)
def pts_of(it):
    if it[0]=='l': return [it[1],it[2]]
    if it[0]=='c': return [it[1],it[4]]
    return []
def objects(pg, r):
    """markers: small closed bezier shapes -> center; lines: long polylines"""
    p=D[pg-1]; M=collections.defaultdict(list); Lns=collections.defaultdict(list)
    for dr in p.get_drawings():
        if not r.intersects(dr["rect"]): continue
        col=dr["color"] or dr["fill"]
        if col is None: continue
        col=tuple(round(c,3) for c in col)
        rc=dr["rect"]; n=len(dr["items"])
        if rc.width<6 and rc.height<6 and all(it[0]=='c' for it in dr["items"]) and n>=4:
            M[col].append(((rc.x0+rc.x1)/2,(rc.y0+rc.y1)/2, rc.width, dr['type']))
        elif n>=8 and 's' in dr['type'] and max(col)>0.05:
            P=[]
            for it in dr["items"]:
                for q in pts_of(it):
                    if not P or abs(P[-1][0]-q.x)+abs(P[-1][1]-q.y)>1e-6: P.append((q.x,q.y))
            Lns[col].append(np.array(P))
    return M, Lns
def mapper(fr, xr, yr, xlog=False, ylog=False):
    x0,x1,y0,y1=fr
    def mx(px):
        f=(np.asarray(px)-x0)/(x1-x0)
        return 10**(np.log10(xr[0])+f*np.log10(xr[1]/xr[0])) if xlog else xr[0]+f*(xr[1]-xr[0])
    def my(py, yrr=yr):
        f=(y1-np.asarray(py))/(y1-y0)
        return 10**(np.log10(yrr[0])+f*np.log10(yrr[1]/yrr[0])) if ylog else yrr[0]+f*(yrr[1]-yrr[0])
    return mx,my
def lin(pa,va,pb,vb,log=False):
    if log:
        la,lb=np.log10(va),np.log10(vb)
        return lambda p: 10**(la+(np.asarray(p)-pa)*(lb-la)/(pb-pa))
    return lambda p: va+(np.asarray(p)-pa)*(vb-va)/(pb-pa)
def assign(markers, lines):
    """index of nearest line (vertical distance at marker x) for each marker"""
    out=[]
    for x,y,*_ in markers:
        dd=[]
        for P in lines:
            o=np.argsort(P[:,0]); dd.append(abs(np.interp(x,P[o,0],P[o,1])-y) if P[:,0].min()-1<=x<=P[:,0].max()+1 else 1e9)
        out.append(int(np.argmin(dd)))
    return np.array(out)
def save(fn, header, cols):
    np.savetxt(fn, np.column_stack(cols), delimiter=",", header=header, comments="", fmt="%.6g")
def ticks_len(pg, r, fr, L=4.0):
    """bottom-axis and left-axis ticks with their lengths"""
    x0,x1,y0,y1=fr; p=D[pg-1]; xt={}; yl={}
    for dr in p.get_drawings():
        if not r.contains(dr["rect"]) or not dr["color"] or max(dr["color"])>0.05: continue
        for it in dr["items"]:
            if it[0]!='l': continue
            a,b=it[1],it[2]
            if abs(a.x-b.x)<0.01 and abs(a.y-b.y)<L and (abs(max(a.y,b.y)-y1)<0.08 or abs(min(a.y,b.y)-y1)<0.08):
                xt[round(a.x,2)]=max(xt.get(round(a.x,2),0),round(abs(a.y-b.y),2))
            if abs(a.y-b.y)<0.01 and abs(a.x-b.x)<L and abs(min(a.x,b.x)-x0)<0.08:
                yl[round(a.y,2)]=max(yl.get(round(a.y,2),0),round(abs(a.x-b.x),2))
    return sorted(xt.items()), sorted(yl.items())
def right_ticks(pg, r, fr, L=4.0):
    x0,x1,y0,y1=fr; p=D[pg-1]; yr={}
    for dr in p.get_drawings():
        if not r.contains(dr["rect"]) or not dr["color"] or max(dr["color"])>0.05: continue
        for it in dr["items"]:
            if it[0]!='l': continue
            a,b=it[1],it[2]
            if abs(a.y-b.y)<0.01 and abs(a.x-b.x)<L and abs(max(a.x,b.x)-x1)<0.08:
                yr[round(a.y,2)]=max(yr.get(round(a.y,2),0),round(abs(a.x-b.x),2))
    return sorted(yr.items())
def black_lines(pg, r, fr, nmin=8):
    x0,x1,y0,y1=fr; p=D[pg-1]; out=[]
    for dr in p.get_drawings():
        if not r.intersects(dr["rect"]) or not dr["color"] or max(dr["color"])>0.05: continue
        if len(dr["items"])<nmin: continue
        P=[]
        for it in dr["items"]:
            for q in pts_of(it):
                if not P or abs(P[-1][0]-q.x)+abs(P[-1][1]-q.y)>1e-6: P.append((q.x,q.y))
        P=np.array(P)
        if P[:,0].max()>x0+2 and P[:,0].min()<x1-2 and P[:,1].min()>y0+1 and P[:,1].max()<y1+1 and np.ptp(P[:,0])>10:
            out.append(P)
    return out


# ---------------------------------------------------------------- figure extraction (Steps 10b, 12a, 12b)
import os
OUT = sys.argv[2] if len(sys.argv) > 2 else "paper_data"
TAG = "extracted 2026-11-18"
SRC = "# read from the PDF's vector drawing (exact plotted positions), {} ({})\n"


def write(name, comment, step, header, arr):
    with open(os.path.join(OUT, name + ".csv"), "w") as o:
        o.write(comment + "\n" + SRC.format(TAG, step) + header + "\n")
        np.savetxt(o, arr, delimiter=",", fmt="%.6g")


def fig9():
    COL = {(0.635, 0.0, 0.0): -3, (0.133, 0.271, 0.408): -2, (0.714, 0.537, 0.004): -1, (0.153, 0.514, 0.294): 0}
    R = fitz.Rect(45, 222, 567, 345)
    PAN = {"a": dict(fr=(66.2, 174.6, 233.6, 312.8), xt=(66.19, 118.3, 170.45), yt=(235.93, 15, 306.87, -45), q="S21_dB",
                     d="Fig. 9(a): Device 2, V_DS = 3 V, |S21| in dB vs frequency"),
           "b": dict(fr=(195.2, 303.9, 232.8, 312.5), xt=(195.18, 247.48, 299.75), yt=(242.74, 30, 302.56, 0), q="gm_mS",
                     d="Fig. 9(b): Device 2, V_DS = 3 V, g_m = Re(Y21) in mS (whole device) vs frequency"),
           "c": dict(fr=(324.6, 433.3, 232.8, 312.5), xt=(324.55, 376.86, 429.12), yt=(251.55, 10, 298.49, 0), q="gds_mS",
                     d="Fig. 9(c): Device 2, V_DS = 3 V, g_ds = Re(Y22) in mS (whole device) vs frequency"),
           "d": dict(fr=(452.6, 562.7, 233.6, 313.5), xt=(452.64, 505.57, 558.5), yt=(236.01, 45, 307.54, -15), q="H21_dB",
                     d="Fig. 9(d): Device 2, V_DS = 3 V, H21 in dB vs frequency")}
    M, Ls = objects(6, R)
    for k, P in PAN.items():
        x0, x1, y0, y1 = P["fr"]
        fx = lin(P["xt"][0], 1.0, P["xt"][2], 100.0, log=True)
        fy = lin(P["yt"][0], P["yt"][1], P["yt"][2], P["yt"][3])
        for col, vg in COL.items():
            mk = np.array(sorted((float(fx(x)), float(fy(y))) for x, y, *_ in M.get(col, []) if x0 - 1 < x < x1 + 1 and y0 - 1 < y < y1 + 1))
            segs = [Pl for Pl in Ls.get(col, []) if Pl[:, 0].max() > x0 + 1 and Pl[:, 0].min() < x1 - 1 and Pl[:, 1].min() > y0 - 2 and Pl[:, 1].max() < y1 + 2]
            Pc = np.vstack(segs); Pc = Pc[(Pc[:, 0] >= x0) & (Pc[:, 0] <= x1)]
            ln = np.column_stack([fx(Pc[:, 0]), fy(Pc[:, 1])]); ln = ln[np.argsort(ln[:, 0])]
            for kind, arr, what in (("measured", mk, "measured (symbols)"), ("paper_model", ln, "the paper's model (line, plotted points)")):
                write(f"fig9{k}_{kind}_{P['q']}_VGS{vg}", f"# Nazir et al., IEEE TED 72(9) 2025, {P['d']}; {what}, V_GS = {vg} V",
                      "Step 10b", f"freq_GHz,{P['q']}", arr)


def fig10_11():
    R = fitz.Rect(45, 357, 567, 482)
    B, RD, G, K = (0.133, 0.271, 0.408), (0.635, 0.0, 0.0), (0.106, 0.337, 0.196), (0.0, 0.0, 0.0)
    PAN = {
        "fig10a": dict(fr=(64.5, 158.7, 378.9, 454.5), x=(68.59, -5, 150.51, 15), xn="Pdel_dBm",
                       s={"Pout_dBm": (B, (454.48, 0, 378.93, 30)), "Gain_dB": (RD, (454.48, 0, 378.93, 30)), "PAE_pct": (K, (454.48, 0, 378.93, 50))},
                       d="Fig. 10(a): Device 1, 10 GHz, V_DS = 15 V, I_DS = 25 mA, Z_S = 25+26j ohm, Z_L = 45+96j ohm; Pout and Gain on the left axis, PAE (%) on the right"),
        "fig10b": dict(fr=(192.8, 283.0, 378.9, 454.5), x=(192.85, -2, 272.96, 6), xn="Pdel_dBm",
                       s={"Pfund_dBm": (RD, (401.59, 0, 386.45, 20)), "IM3_dBm": (B, (401.59, 0, 386.45, 20)), "OIP3_over_PDC_dB": (G, (397.8, 0, 378.89, 20))},
                       d="Fig. 10(b): Device 1, two-tone, f0 = 10 GHz, delta f = 10 MHz, V_DS = 15 V, Z_S = 20+26j ohm, Z_L = 50+140j ohm; P_FUND and IM3 left axis, OIP3/P_DC right axis"),
        "fig11a": dict(fr=(326.3, 435.6, 378.9, 455.7), x=(326.31, 1, 386.64, 10), xlog=True, xn="freq_GHz",
                       s={"H21_dB": (RD, (455.7, 0, 378.9, 45)), "U_dB": (B, (455.7, 0, 378.9, 45))},
                       d="Fig. 11(a): Device 1, V_GS = 0 V, V_DS = 10 V, current gain H21 and unilateral power gain U in dB vs frequency"),
        "fig11b": dict(fr=(454.4, 562.6, 378.9, 454.5), x=(454.36, 0, 559.07, 600), xn="IDS_mA_per_mm",
                       s={"fT_GHz": (RD, (454.48, 0, 378.93, 60)), "fMAX_GHz": (B, (454.48, 0, 378.93, 60))},
                       d="Fig. 11(b): Device 1, V_DS = 10 V, f_T and f_MAX vs I_DS (mA/mm)")}
    M, Ls = objects(6, R)
    for k, P in PAN.items():
        x0, x1, y0, y1 = P["fr"]
        fx = lin(P["x"][0], P["x"][1], P["x"][2], P["x"][3], log=P.get("xlog", False))
        for q, (col, yt) in P["s"].items():
            fy = lin(*yt)
            mk = np.array(sorted((float(fx(x)), float(fy(y))) for x, y, w, t in M.get(col, []) if x0 - 1 < x < x1 + 1 and y0 - 1 < y < y1 + 1 and 2.5 < w < 2.9))
            segs = black_lines(6, R, P["fr"]) if col == K else [Pl for Pl in Ls.get(col, []) if Pl[:, 0].max() > x0 + 1 and Pl[:, 0].min() < x1 - 1 and Pl[:, 1].min() > y0 - 2 and Pl[:, 1].max() < y1 + 9]
            Pc = np.vstack(segs); Pc = Pc[(Pc[:, 0] >= x0) & (Pc[:, 0] <= x1) & (Pc[:, 1] <= y1 + 0.01)]
            ln = np.column_stack([fx(Pc[:, 0]), fy(Pc[:, 1])]); ln = ln[np.argsort(ln[:, 0])]
            for kind, arr, what in (("measured", mk, "measured (symbols)"), ("paper_model", ln, "the paper's model (line, plotted points)")):
                write(f"{k}_{kind}_{q}", f"# Nazir et al., IEEE TED 72(9) 2025, {P['d']}; {what}", "Step 10b", f"{P['xn']},{q}", arr)


def fig7c():
    R = fitz.Rect(390, 345, 567, 488)
    M, Ls = objects(5, R)
    TH = {(0.635, 0.0, 0.0): 0.1, (0.133, 0.271, 0.408): 0.2, (0.714, 0.537, 0.004): 0.3, (0.153, 0.514, 0.294): 0.4, (0.271, 0.271, 0.404): 0.5}
    fx = lin(417.01, -3, 535.23, 1); fgm = lin(457.3, 0, 354.14, 500); fvs = lin(459.38, 0, 354.14, 500)
    for col, th in TH.items():
        a, b = Ls[col]
        gm, vs = (a, b) if a[:, 1].max() > b[:, 1].max() else (b, a)
        for P, fy, q, unit in ((gm, fgm, "GM", "mS"), (vs, fvs, "VSAT", "km_per_s")):
            d = np.column_stack([fx(P[:, 0]), fy(P[:, 1])]); d = d[np.argsort(d[:, 0])]
            write(f"fig7c_paper_model_{q}_th2_{th:.1f}",
                  f"# Nazir et al., IEEE TED 72(9) 2025, Fig. 7(c): simulation study, theta_SAT2 = {th}, theta_SAT1 = 0.3, V_SAT0 = 1e5 m/s, n_max = 5e16 m^-2; the paper's model line ({q})",
                  "Step 12a", f"VGS_V,{q}_{unit}", d)


def fig8():
    R = fitz.Rect(45, 40, 567, 196)
    M, Ls = objects(6, R)
    COL = {(0.635, 0.0, 0.0): -3, (0.133, 0.271, 0.408): -2, (0.714, 0.537, 0.004): -1, (0.153, 0.514, 0.294): 0}
    # chart centre and rim from the chart outlines (Smith: unit circle; polar: outer circle = rim value)
    PAN = {"S11": (105.95, 123.45, 44.95, 1.0, (55, 156)), "S12": (245.05, 126.65, 49.95, 0.6, (190, 300)),
           "S21": (385.35, 126.25, 49.95, 3.0, (330, 440)), "S22": (504.65, 121.15, 44.95, 1.0, (455, 555))}
    for q, (cx, cy, Rr, sc, (xa, xb)) in PAN.items():
        for col, vg in COL.items():
            mk = np.array([(x, y) for x, y, w, t in M[col] if xa < x < xb and 70 < y < 180 and 2.5 < w < 4.5])
            ln = [P for P in Ls[col] if xa < P[:, 0].mean() < xb][0]
            chart = "Smith chart, |Gamma| = 1 at the rim" if q in ("S11", "S22") else f"polar chart, rim = {sc}"
            for kind, arr, what in (("measured", mk, "measured (symbols, 1-110 GHz, order along the trace)"), ("paper_model", ln, "the paper's model (line)")):
                g = ((arr[:, 0] - cx) - 1j * (arr[:, 1] - cy)) / Rr * sc
                write(f"fig8_{kind}_{q}_VGS{vg}", f"# Nazir et al., IEEE TED 72(9) 2025, Fig. 8: Device 2, V_DS = 3 V, {q}, V_GS = {vg} V; {what}; {chart}",
                      "Step 12b", "Re,Im", np.column_stack([g.real, g.imag]))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    fig7c(); fig8(); fig9(); fig10_11()
    print("wrote", len([f for f in os.listdir(OUT) if f.startswith(("fig7c", "fig8_", "fig9", "fig10", "fig11"))]), "CSVs to", OUT)
