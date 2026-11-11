"""Rebuild exp3_{TP}_peakfeats.npz caches from corrected exp3_load labels (canonical split-Gaussian
12-(n,m) fits per well). Run after any plate-label change."""
import os, sys, warnings, importlib.util
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[v]="1"
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from repo_paths import data_path, results_dir, env_path
import numpy as np, exp3_load as L
PD=os.path.join(os.path.dirname(os.path.abspath(__file__)),"vendored")
# rebuilt caches are written under results/, never over the shipped cache in data/
OUT=env_path("EXP3_OUT", "") or results_dir("04_drug_response_prediction", "rebuilt_cache"); os.makedirs(OUT, exist_ok=True)
spec=importlib.util.spec_from_file_location("cie",os.path.join(PD,"chirality_interpolated_extraction.py"))
cie=importlib.util.module_from_spec(spec); spec.loader.exec_module(cie)
CHIR=cie.CHIRALITY_POSITIONS; EXC=np.array([570,640,670,750,780],float); exmin,exmax=EXC.min(),EXC.max()
GROUPS=["sham","tbi","resv_tbi","mg132_tbi","resv_sham","mg132_sham","swcnt_noncm"]
PLATE={"sham":"SHAM","resv_sham":"SHAM","mg132_sham":"SHAM","swcnt_noncm":"SHAM",
       "tbi":"TBI","resv_tbi":"TBI","mg132_tbi":"TBI"}
def well_chir(mat,em):
    sp,emf,exf=cie.build_interpolated_surface(mat,em,EXC,factor=cie.INTERP_FACTOR)
    mx=np.zeros(12); au=np.zeros(12)
    for k,(name,c) in enumerate(CHIR.items()):
        tex=min(max(c["excitation_nm"],exmin),exmax)
        be,bx,m=cie.find_local_max_in_radius(sp,emf,exf,c["emission_nm"],tex,cie.SEARCH_RADIUS_NM)
        ep,ip=cie.extract_consistent_gradient_points(sp,emf,be,bx)
        d=cie.fit_gaussian_and_extract(ep,ip,be,m); mx[k]=max(d["gauss_max"],0); au[k]=max(np.nan_to_num(d["gauss_auc"]),0)
    return mx,au
for tp in ("3H","6H"):
    cs,em=L.load_plate("SHAM",tp); ct,_=L.load_plate("TBI",tp); Gs=L.groups("SHAM"); Gt=L.groups("TBI")
    out={}
    for g in GROUPS:
        cube=cs if PLATE[g]=="SHAM" else ct; idx=(Gs if PLATE[g]=="SHAM" else Gt)[g]
        MX=np.zeros((len(idx),12)); AU=np.zeros((len(idx),12))
        for i,w in enumerate(idx): MX[i],AU[i]=well_chir(cube[w].astype(float),em)
        out[f"MX_{g}"]=MX; out[f"AU_{g}"]=AU
    p=os.path.join(OUT,"exp3_%s_peakfeats.npz"%tp); np.savez(p,**out); print("rebuilt",p,flush=True)
print("done")
