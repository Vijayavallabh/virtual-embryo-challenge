# RECONSTRUCTED: lines 2-25 are the exact body of the quoted heredoc that wrote /tmp/lvl.py (deleted after the run); run once as: cd /mnt/md0/IITM/BackUp/Home/vijayavallabh/virtual-embryo-challenge && timeout 420 .venv/bin/python -I /tmp/lvl.py 2>&1 | tail -8
import sys, numpy as np, anndata as ad
from scipy.stats import spearmanr
def prop_mean(p):
    a = ad.read_h5ad(p)
    X = np.asarray(a.X.todense() if hasattr(a.X,'todense') else a.X, dtype=np.float64)
    P = np.expm1(X)/1e4
    return P.mean(0), [str(g) for g in a.var_names], a.n_obs
wt, gw, nw = prop_mean('data/raw/T2/E8.75.h5ad')
print('WT E8.75 cells', nw, 'genes', len(gw))
pb_ref = np.log1p(1e4*wt).astype(float)   # reference pseudobulk proxy, gene level
for name, path in [('M','data/raw/T3/E9.5_mab21l2_ko.h5ad'), ('G1','data/raw/T3/E8.75_gata4_ko_1.h5ad')]:
    ko, gk, nk = prop_mean(path)
    assert gk == gw, 'gene order differs'
    lr = np.log((ko+1e-5)/(wt+1e-5))
    lr = lr - np.average(lr, weights=wt)          # den-weighted total held, as gene_factor does
    r = spearmanr(lr, pb_ref).statistic
    ra = spearmanr(np.abs(lr), pb_ref).statistic
    print(name, 'cells', nk, '| SD(ln phi) %.3f' % lr.std(),
          '| rho(ln phi, pb_ref) %+.3f' % r, '| rho(|ln phi|, pb_ref) %+.3f' % ra)
    np.save('/tmp/lr_%s.npy' % name, lr)
m = np.load('/tmp/lr_M.npy'); g = np.load('/tmp/lr_G1.npy')
print('rho(ln phi_M, ln phi_G1) %+.3f | rho(|M|,|G1|) %+.3f' % (spearmanr(m,g).statistic, spearmanr(np.abs(m),np.abs(g)).statistic))
c2 = 0.82*(m+g)
print('C2-like: rho(ln phi_C2, pb_ref) %+.3f' % spearmanr(c2, pb_ref).statistic)
