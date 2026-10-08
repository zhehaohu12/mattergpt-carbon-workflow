#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare five groups of generated pure-carbon structures with a C32-C2/m
reference using the Oganov–Valle fingerprint cosine distance.

Reference:
A. R. Oganov and M. Valle, J. Chem. Phys. 130, 104504 (2009).
DOI: 10.1063/1.3079326

Pure-carbon fingerprint: F(R)=g(R)-1
Distance: D=0.5*[1-cos(Fref,Fi)]
Similarity for plotting: S=1-D

Default parameters: Rmax=15 Å, dR=0.05 Å, sigma=0.075 Å

Expected files:
resultsC0000.csv ... resultsC9000.csv
Default reference: C32-C2-m-no-opt.vasp

Groups:
G1=C0000+C1000
G2=C2000+C3000
G3=C4000+C5000
G4=C6000+C7000
G5=C8000+C9000

Usage:
python compare_five_groups_to_C32_oganov.py .
python compare_five_groups_to_C32_oganov.py . --reference your_C32.vasp
"""

import os
os.environ.setdefault("OMP_NUM_THREADS","1")
os.environ.setdefault("MKL_NUM_THREADS","1")
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")

import argparse, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from pymatgen.io.vasp.inputs import Poscar

RUNS=list(range(0,10000,1000))
GROUPS=[
    ("Group 1",[0,1000]),
    ("Group 2",[2000,3000]),
    ("Group 3",[4000,5000]),
    ("Group 4",[6000,7000]),
    ("Group 5",[8000,9000]),
]

def find_files(base):
    found={}
    for p in base.glob("resultsC*.csv"):
        m=re.search(r"resultsC(\d{4})",p.name,re.I)
        if m: found[int(m.group(1))]=p
    miss=[x for x in RUNS if x not in found]
    if miss:
        raise FileNotFoundError("Missing: "+", ".join(f"resultsC{x:04d}.csv" for x in miss))
    return found

def detect_col(df, requested=None):
    df.columns=[str(c).replace("\ufeff","").strip() for c in df.columns]
    mp={c.lower():c for c in df.columns}
    if requested and requested.lower() in mp: return mp[requested.lower()]
    for k in ["poscar","structure","vasp","contcar","poscar_text","generated_poscar","structure_poscar"]:
        if k in mp: return mp[k]
    raise ValueError(f"Cannot detect POSCAR column. Columns={list(df.columns)}")

def parse_struct(v):
    s=str(v)
    if "\\n" in s: s=s.replace("\\n","\n")
    return Poscar.from_str(s).structure

def pure_c(s):
    return {e.symbol for e in s.composition.elements}=={"C"}

def fingerprint(s,rmax=15.0,dr=0.05,sigma=0.075):
    N=len(s); V=float(s.volume)
    if N<=0 or V<=0: raise ValueError("Invalid structure")
    r=np.arange(0.0,rmax+0.5*dr,dr,dtype=float)
    ds=[]
    for site in s:
        for n in s.get_neighbors(site,rmax):
            d=float(n.nn_distance)
            if d>1e-10: ds.append(d)
    ds=np.asarray(ds,float)
    if ds.size==0: raise ValueError("No pairs within Rmax")
    g=np.zeros_like(r)
    norm_pair=(N**2)/V
    for st in range(0,len(ds),2000):
        d=ds[st:st+2000]
        diff=r[:,None]-d[None,:]
        gauss=np.exp(-0.5*(diff/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
        w=1.0/(4*np.pi*(d**2)*norm_pair)
        g += gauss @ w
    F=g-1.0
    nrm=np.linalg.norm(F)
    if nrm<=1e-15: raise ValueError("Zero-norm fingerprint")
    return (F/nrm).astype(np.float32)

def compare(ref,fp):
    c=float(np.clip(np.dot(ref,fp),-1,1))
    D=0.5*(1-c)
    return c,D,1-D

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("directory",nargs="?",default=".")
    ap.add_argument("--reference",default="C32-C2-m-no-opt.vasp")
    ap.add_argument("--poscar-column",default=None)
    ap.add_argument("--rmax",type=float,default=15.0)
    ap.add_argument("--dr",type=float,default=0.05)
    ap.add_argument("--sigma",type=float,default=0.075)
    ap.add_argument("--bin-width",type=float,default=0.025)
    a=ap.parse_args()

    base=Path(a.directory).resolve()
    refp=base/a.reference
    if not refp.exists(): raise FileNotFoundError(refp)
    files=find_files(base)

    refs=Poscar.from_file(str(refp)).structure
    if not pure_c(refs): raise ValueError("Reference is not pure carbon")
    print(f"Reference: {refp.name}",flush=True)
    print("Calculating reference fingerprint...",flush=True)
    fref=fingerprint(refs,a.rmax,a.dr,a.sigma)

    edges=np.arange(0,1+a.bin_width,a.bin_width,dtype=float)
    edges[-1]=1.0
    centers=(edges[:-1]+edges[1:])/2
    ddf=pd.DataFrame({"distance_left":edges[:-1],"distance_right":edges[1:],"distance_bin_center":centers})
    sdf=pd.DataFrame({"similarity_left":edges[:-1],"similarity_right":edges[1:],"similarity_bin_center":centers})

    allrows=[]; sums=[]; invalid=[]; probs=[]

    for gname,ids in GROUPS:
        Ds=[]; Ss=[]
        print(f"\n{gname}",flush=True)
        for rid in ids:
            cp=files[rid]
            df=pd.read_csv(cp)
            col=detect_col(df,a.poscar_column)
            good=0
            for idx,v in df[col].items():
                try:
                    s=parse_struct(v)
                    if not pure_c(s):
                        invalid.append([gname,cp.name,int(idx),"not pure carbon"]); continue
                    fp=fingerprint(s,a.rmax,a.dr,a.sigma)
                    c,D,S=compare(fref,fp)
                    Ds.append(D); Ss.append(S)
                    allrows.append([gname,cp.name,int(idx),len(s),c,D,S])
                    good+=1
                except Exception as e:
                    invalid.append([gname,cp.name,int(idx),str(e)])
            print(f"  {cp.name}: valid pure-C={good}",flush=True)

        Ds=np.asarray(Ds,float); Ss=np.asarray(Ss,float)
        if Ds.size==0: raise RuntimeError(f"No valid structures in {gname}")

        dc,_=np.histogram(Ds,bins=edges); sc,_=np.histogram(Ss,bins=edges)
        dp=dc/dc.sum(); sp=sc/sc.sum()
        ddf[gname.replace(" ","_")+"_pct"]=dp*100
        sdf[gname.replace(" ","_")+"_pct"]=sp*100
        probs.append(dp)

        sums.append({
            "group":gname,
            "n_valid_structures":len(Ds),
            "mean_D_to_C32":Ds.mean(),
            "median_D_to_C32":np.median(Ds),
            "std_D_to_C32":Ds.std(ddof=1),
            "min_D_to_C32":Ds.min(),
            "max_D_to_C32":Ds.max(),
            "mean_S_to_C32":Ss.mean(),
            "median_S_to_C32":np.median(Ss),
            "std_S_to_C32":Ss.std(ddof=1),
            "fraction_S_ge_0.90_percent":(Ss>=0.90).mean()*100,
            "fraction_S_ge_0.95_percent":(Ss>=0.95).mean()*100,
        })
        print(f"  -> n={len(Ds)}, mean D={Ds.mean():.6f}, mean S={Ss.mean():.6f}",flush=True)

    gcols=[f"Group_{i}_pct" for i in range(1,6)]
    for df in [ddf,sdf]:
        df["five_group_mean_pct"]=df[gcols].mean(axis=1)
        df["five_group_sd_pct"]=df[gcols].std(axis=1,ddof=1)

    summary=pd.DataFrame(sums)
    js=[]; jsrows=[]
    for i in range(5):
        for j in range(i+1,5):
            v=float(jensenshannon(probs[i],probs[j],base=2.0)**2)
            js.append(v)
            jsrows.append([f"Group {i+1}",f"Group {j+1}",v])

    mD=summary["mean_D_to_C32"].mean(); sdD=summary["mean_D_to_C32"].std(ddof=1)
    mS=summary["mean_S_to_C32"].mean(); sdS=summary["mean_S_to_C32"].std(ddof=1)
    overall=pd.DataFrame([{
        "reference_file":refp.name,
        "mean_of_group_mean_D_to_C32":mD,
        "SD_of_group_mean_D_to_C32":sdD,
        "CV_of_group_mean_D_percent":sdD/mD*100 if mD else np.nan,
        "mean_of_group_mean_S_to_C32":mS,
        "SD_of_group_mean_S_to_C32":sdS,
        "CV_of_group_mean_S_percent":sdS/mS*100 if mS else np.nan,
        "mean_pairwise_JS_divergence":np.mean(js),
        "max_pairwise_JS_divergence":np.max(js),
    }])

    pd.DataFrame(allrows,columns=["group","source_file","source_row","n_atoms","oganov_cosine_to_C32","oganov_distance_to_C32_D","similarity_to_C32_1_minus_D"]).to_csv(base/"01_all_structures_to_C32_oganov.csv",index=False,float_format="%.8f")
    summary.to_csv(base/"02_five_group_C32_summary.csv",index=False,float_format="%.8f")
    ddf.to_csv(base/"03_C32_distance_distribution_five_groups.csv",index=False,float_format="%.8f")
    sdf.to_csv(base/"04_C32_similarity_distribution_five_groups.csv",index=False,float_format="%.8f")
    pd.DataFrame(jsrows,columns=["group_i","group_j","JS_divergence_C32_distance_distribution"]).to_csv(base/"05_C32_five_group_JS_divergence.csv",index=False,float_format="%.8f")
    overall.to_csv(base/"06_C32_five_group_overall.csv",index=False,float_format="%.8f")
    if invalid:
        pd.DataFrame(invalid,columns=["group","source_file","row","reason"]).to_csv(base/"07_C32_invalid_rows.csv",index=False)

    for mode,df,xcol,xlab,prefix in [
        ("D",ddf,"distance_bin_center","Oganov fingerprint cosine distance to C32-C2/m, D","08_C32_distance_distribution_five_groups"),
        ("S",sdf,"similarity_bin_center","Structural fingerprint similarity to C32-C2/m, S = 1 - D","09_C32_similarity_distribution_five_groups"),
    ]:
        fig,ax=plt.subplots(figsize=(8.5,5.7))
        for i in range(1,6):
            ax.plot(df[xcol],df[f"Group_{i}_pct"],marker="o",markersize=3,linewidth=1.5,label=f"Group {i}")
        ax.set_xlabel(xlab)
        ax.set_ylabel("Fraction of generated structures (%)")
        ax.set_xlim(0,1); ax.legend(); ax.grid(alpha=0.2); fig.tight_layout()
        fig.savefig(base/f"{prefix}.png",dpi=300,bbox_inches="tight")
        fig.savefig(base/f"{prefix}.pdf",bbox_inches="tight")
        plt.close(fig)

    (base/"10_C32_oganov_method_settings.txt").write_text(
        f"""Reference: Oganov & Valle, J. Chem. Phys. 130, 104504 (2009)
Reference structure: {refp.name}
Fingerprint: F(R)=g(R)-1
Distance: D=0.5*[1-cos(Fref,Fi)]
Similarity for plotting: S=1-D
Rmax={a.rmax} Å
dR={a.dr} Å
sigma={a.sigma} Å
Each group's histogram is normalized by its own number of valid structures.
This is a continuous fingerprint-similarity analysis, not a strict equivalence test.
""",encoding="utf-8")

    print("\nGROUP SUMMARY",flush=True)
    print(summary.to_string(index=False,float_format=lambda x:f"{x:.6f}"),flush=True)
    print("\nOVERALL",flush=True)
    print(overall.to_string(index=False,float_format=lambda x:f"{x:.6f}"),flush=True)
    print("\nDone.",flush=True)

if __name__=="__main__":
    main()
