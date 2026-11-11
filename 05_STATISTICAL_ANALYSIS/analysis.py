"""Spatial, anomaly, variability, state, Mann–Kendall and Sen analyses."""
from pathlib import Path
import yaml, json
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import geometry_mask
from scipy.stats import norm, theilslopes

ROOT=Path(__file__).resolve().parents[1]; C=yaml.safe_load((ROOT/'00_CONFIGURATION/config.yaml').read_text())['study']
OUT=ROOT/'08_OUTPUTS'; R=OUT/'Rasters'; T=OUT/'Tables'; T.mkdir(parents=True,exist_ok=True)
def mk(series):
    x=np.asarray(series,dtype=float); n=x.size; s=0
    for k in range(n-1): s+=np.sign(x[k+1:]-x[k]).sum()
    _,counts=np.unique(x,return_counts=True); tie=np.sum(counts*(counts-1)*(2*counts+5))
    var=(n*(n-1)*(2*n+5)-tie)/18
    if var<=0:return 0.,1.
    z=(s-np.sign(s))/np.sqrt(var) if s else 0.
    return float(z),float(2*norm.sf(abs(z)))
def run(data=None,years=None,profile=None):
    if data is None:
        q=np.load(R/'january_stack.npz'); data=q['data']; years=q['years'].tolist()
        with rasterio.open(R/f"jan_{years[0]}_c.tif") as d: profile=d.profile
    years=np.asarray(years); pyears=(years>=C['primary_start_year'])&(years<=C['primary_end_year']); cyears=(years>=C['climatology_start_year'])&(years<=C['climatology_end_year'])
    y=years[pyears]; x=data[pyears].astype('float64'); nodata=C['nodata']; x[x==nodata]=np.nan
    base=np.nanmean(np.where(np.isfinite(data[cyears])&(data[cyears]!=nodata),data[cyears],np.nan),axis=0)
    anom=x-base; mean=np.nanmean(x,axis=0); sd=np.nanstd(x,axis=0,ddof=1); cv=sd/np.abs(mean)*100
    h,w=mean.shape; slope=np.full((h,w),np.nan); zmap=slope.copy(); pmap=slope.copy()
    for i in range(h):
        for j in range(w):
            s=x[:,i,j]; ok=np.isfinite(s)
            if ok.sum()<3: continue
            z,p=mk(s[ok]); zmap[i,j]=z;pmap[i,j]=p
            slope[i,j]=theilslopes(s[ok],y[ok])[0]*10
    sig=C['significance_alpha']; cat=np.where((pmap<sig)&(slope>0),1,np.where((pmap<sig)&(slope<0),-1,0)).astype('int8')
    prof=profile.copy(); prof.update(count=1,compress='deflate')
    def save(name,a,dtype='float32',nd=nodata):
        pp=prof.copy();pp.update(dtype=dtype,nodata=nd)
        with rasterio.open(R/f'{name}.tif','w',**pp) as d:d.write(np.where(np.isfinite(a),a,nd).astype(dtype),1)
    save('jan_period_mean_c',mean); save('jan_period_sd_c',sd); save('jan_period_cv_pct',cv); save('jan_climatology_configured_baseline_c',base)
    save('sen_slope_c_per_decade',slope);save('mann_kendall_z',zmap);save('mann_kendall_p',pmap);save('trend_class',cat,'int8',-128)
    for i,yr in enumerate(y):save(f'anomaly_{yr}_c',anom[i])
    state=gpd.read_file(ROOT/C['boundary']); field=C['state_field']; state=state[state[field].isin(C['states'])].to_crs(prof['crs'])
    rows=[]; areas=[]; transform=prof['transform']
    # Regional and state annual means use equal-area weighting via cell-area weights.
    from rasterio.transform import xy
    lat=np.array([xy(transform,0,j)[1] for j in range(h)])[:,None]
    weights=np.cos(np.deg2rad(lat)); valid=np.isfinite(x)
    reg=np.nansum(np.where(valid,x*weights,np.nan),axis=(1,2))/np.nansum(np.where(valid,weights,np.nan),axis=(1,2))
    for yr,v in zip(y,reg):rows.append({'year':int(yr),'geography':'Southwestern Nigeria','mean_c':v})
    pixarea=abs(transform.a*transform.e)
    for _,f in state.iterrows():
        mask=geometry_mask([f.geometry.__geo_interface__],(h,w),transform,invert=True)
        for k,yr in enumerate(y):
            a=x[k]; vals=a[mask & np.isfinite(a)]; rows.append({'year':int(yr),'geography':str(f[field]),'mean_c':float(np.mean(vals)) if vals.size else np.nan})
        sv=slope[mask & np.isfinite(slope)]; pv=pmap[mask & np.isfinite(pmap)]; cvv=cat[mask]
        areas.append({'state':str(f[field]),'sen_slope_mean_c_decade':float(np.mean(sv)) if sv.size else np.nan,'sen_slope_median_c_decade':float(np.median(sv)) if sv.size else np.nan,'significant_warming_pct':float(np.mean(cvv==1)*100),'significant_cooling_pct':float(np.mean(cvv==-1)*100),'not_significant_pct':float(np.mean(cvv==0)*100)})
    ts=pd.DataFrame(rows);ts.to_csv(T/'regional_state_annual_temperature.csv',index=False)
    pd.DataFrame(areas).to_csv(T/'state_trend_statistics.csv',index=False)
    summary=ts.groupby('geography').mean_c.agg(['mean','std','min','max']).reset_index();summary.to_csv(T/'temperature_descriptive_statistics.csv',index=False)
    pd.DataFrame({'year':y,'regional_mean_c':reg,'regional_anomaly_c':reg-np.nanmean(base)}).to_csv(T/'regional_anomalies.csv',index=False)
    counts=pd.Series(cat[cat!=-128] if False else cat.ravel()).value_counts();pd.DataFrame([{'class':k,'pixel_share_pct':float(counts.get(k,0)/counts.sum()*100)} for k in [-1,0,1]]).to_csv(T/'trend_area_share.csv',index=False)
    print('Analysis outputs written to',OUT)
if __name__=='__main__':run()
